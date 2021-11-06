"""
    train_prediction_model/__main__.py
    Dylan H. Ross

    this script coordinates training of a CCS prediction model
 Trains a CCS prediction model using K-means clustering and individual SVR models for each cluster
"""

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVR
import pickle
import json
from matplotlib import rcParams
from numpy import mean, median, abs, sum, cumsum, histogram, sqrt
import sys
from sqlite3 import connect

from .C3SData.data import C3SD
from .kmcm import KMCMulti, kmcm_p_grid
from .config import seed, gs_n_jobs, n_clusters, C, gamma, hyperparam_permute


# set the font size for plots
rcParams['font.size'] = 6


def metrics(model, data):
    """
metrics
    description:
        computes a standard set of performance metrics using a trained model and dataset
            * training and test set R-squared (R2)
            * training and test set root mean squared error (RMSE)
            * training and test set mean absolute error (MAE)
            * cumulative error distribution at the <1, <3, <5, and <10% levels
              for training and test set (CE135A)
"""
    summary = {}
    for s in ['train', 'test']:
        if s == 'train':
            y = data.y_train_
            y_pred = model.predict(data.X_train_ss_)
        if s == 'test':
            y = data.y_test_
            y_pred = model.predict(data.X_test_ss_)

        # compute metrics
        abs_y_err = abs(y_pred - y)
        r2 = r2_score(y, y_pred)
        mae = mean(abs_y_err)
        mdae = median(abs_y_err)
        rmse = sqrt(mean_squared_error(y, y_pred))
        y_err_percent = 100. * abs_y_err / y
        cum_err = cumsum(histogram(y_err_percent, [_ for _ in range(101)])[0])
        cum_err = 100. * cum_err / sum(cum_err)
        ce1, ce3, ce5, ceA = cum_err[0], cum_err[2], cum_err[4], cum_err[9]
        summary[s] = {'R2': r2, 'MAE': mae, 'MDAE': mdae, 'RMSE': rmse, 'CE135A': [ce1, ce3, ce5, ceA]}
    return summary


def summary_figure(summary, r2_range=[0.95, 1.]):
    """
summary_figure
    description:
        produces a summary figure displaying the results from a trial
    parameters:
        [r2_range (list(float))] -- upper and lower bounds of R-squared y axis [optional, default=[0.95, 1.]]
        [save (None or str)] -- if a filename is provided, save the figure (optional, default=None)
"""
    fig = plt.figure(figsize=(3.33, 1.75))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 3, 5])

    # R-squared
    ax1 = fig.add_subplot(gs[0])
    rsq_trn = summary['train']['R2']
    rsq_tst = summary['test']['R2']
    w1 = 0.15
    ax1.bar([1 - w1 / 2, 1 + w1 / 2], [rsq_trn, rsq_tst], color=['b', 'r'], width=w1)
    for d in ['top', 'right']:
        ax1.spines[d].set_visible(False)
    ax1.set_xticks([])
    ax1.set_ylabel(r'R$^2$')
    ax1.set_ylim(r2_range)
    ax1.set_xlim([0.75, 1.25])

    # MAE, MDAE and RMSE
    ax2 = fig.add_subplot(gs[1])
    mae_trn = summary['train']['MAE']
    mae_tst = summary['test']['MAE']
    mdae_trn = summary['train']['MDAE']
    mdae_tst = summary['test']['MDAE']
    mse_trn = summary['train']['RMSE']
    mse_tst = summary['test']['RMSE']
    ax2.bar([0.875, 1.125], [mae_trn, mae_tst], color=['b', 'r'], width=0.25)
    ax2.bar([1.875, 2.125], [mdae_trn, mdae_tst], color=['b', 'r'], width=0.25)
    ax2.bar([2.875, 3.125], [mse_trn, mse_tst], color=['b', 'r'], width=0.25)
    for d in ['top', 'right']:
        ax2.spines[d].set_visible(False)
    ax2.set_xticks([1, 2, 3])
    ax2.set_xticklabels(['MAE', 'MDAE', 'RMSE'], rotation='vertical')
    ax2.set_ylabel(r'CCS (Å$^2$)')

    # CE135A
    ax3 = fig.add_subplot(gs[2])
    x1 = [_ - 0.125 for _ in range(1, 5)]
    y1 = [100. * summary['train']['CE135A'][i] for i in range(4)]
    x2 = [_ + 0.125 for _ in range(1, 5)]
    y2 = [100. * summary['test']['CE135A'][i] for i in range(4)]

    ax3.bar(x1, y1, color='b', width=0.25)
    ax3.bar(x2, y2, color='r', width=0.25)

    for d in ['top', 'right']:
        ax3.spines[d].set_visible(False)
    ax3.set_xlabel('pred. error (%)')
    ax3.set_xticks([1, 2, 3, 4])
    ax3.set_xticklabels(['<1', '<3', '<5', '<10'])
    ax3.set_ylabel('proportion (%)')

    plt.tight_layout()
    plt.savefig('kmcm_svr_final_metrics.png', dpi=400, bbox_inches='tight')


def add_predicted_ccs(model, data, cur):
    """ adds predicted CCS to the predicted table in C3S.db """
    # clear out the predicted table if there are already values in it
    cur.execute("""DELETE FROM predicted""")

    # retrieve attributes that are already stored in the C3SD instance
    gids, names, mzs, adducts, smis = data.g_id_, data.cmpd_, data.mz_, data.adduct_, data.smi_

    # generate predicted CCS, class labels, and prediction errors for each compound
    Xss = data.SScaler_.transform(data.X_)
    pred_ccss = model.predict(Xss)
    pred_errors = pred_ccss - data.y_
    class_labels = model.kmeans_.predict(Xss)

    qry = """INSERT INTO predicted VALUES (?,?,?,?,?,?,?,?,?)"""
    for gid, name, adduct, mz, pred_ccs, smi, class_label, pred_error in zip(gids, names, adducts, mzs, pred_ccss, 
                                                                             smis, class_labels, pred_errors):
        qdata = (gid, name, adduct, mz, pred_ccs, smi, int(class_label), pred_error, None)
        cur.execute(qry, qdata)


def main(db_path):
    data = C3SD(db_path, seed=seed)
    data.featurize()
    data.train_test_split('ccs')
    data.center_and_scale()
    kmcm_svr_p_grid = kmcm_p_grid(n_clusters, {'C': C, 'gamma': gamma}, permute=hyperparam_permute)

    kmcm_svr_gs = GridSearchCV(KMCMulti(seed=seed, use_estimator=SVR(cache_size=1024, tol=5e-4)),
                               param_grid=kmcm_svr_p_grid, n_jobs=gs_n_jobs, cv=5, scoring='neg_mean_squared_error', verbose=3)
    kmcm_svr_gs.fit(data.X_train_ss_, data.y_train_)
    kmcm_svr_best = kmcm_svr_gs.best_estimator_
    print(kmcm_svr_gs.best_params_)

    summary = metrics(kmcm_svr_best, data)

    # dump a summary of performance metrics to a .json file
    with open('kmcm_svr_final_metrics.json', 'w') as f:
        json.dump(summary, f, indent=4)

    summary_figure(summary)

    with open('kmcm_svr_final.pickle', 'wb') as pf:
        pickle.dump(kmcm_svr_best, pf)

    with open('OHEncoder.pickle', 'wb') as ohe, open('LEncoder.pickle', 'wb') as le, open('SScaler.pickle', 'wb') as ss:
        pickle.dump(data.OHEncoder_, ohe)
        pickle.dump(data.LEncoder_, le)
        pickle.dump(data.SScaler_, ss)

    # add predicted CCS values to predicted table in database
    con = connect('C3S.db')
    cur = con.cursor()
    add_predicted_ccs(kmcm_svr_best, data, cur)
    con.commit()
    con.close()


if __name__ == '__main__':
    """ main execution sequence """

    # make sure that a path to C3S.db is provided
    if len(sys.argv) < 2:
        msg = "train_prediction_model: __main__: path to C3S.db must be provided"
        raise ValueError(msg)

    # train the CCS prediction model using the provided CCS database
    main(sys.argv[1])

