"""
    generate_mqns.py
    Dylan Ross
    2018/10/19

        Utility to generate MQNs as molecular descriptors
"""

from sqlite3 import connect
from rdkit import Chem
from rdkit.Chem import Descriptors


def get_MQNs(smi):
    """
get_MQNs
    description:
        Computes the complete set of 42 MQNs as described in:
            Nguyen et al. ChemMedChem 4:1803-5 (2009)
    parameters:
        smi (str) -- SMILES string of structure
    returns:
        (list(float) or None) -- array of features, None if anything goes wrong
"""
    try:
        features = Descriptors.rdMolDescriptors.MQNs_(Chem.MolFromSmiles(smi))
        return features
    except Exception as e:
        #print(e)
        print("\tunable to generate MQNs for SMILES: '{}'".format(smi))
        return None


def main():
    """ main execution sequence """

    con = connect("C3S.db")
    cur = con.cursor()

    print("Computing MQNs for all compounds with SMILES structures ...")

    # generate the rdk features and MQNs
    qry = "SELECT g_id, smi FROM master WHERE smi IS NOT NULL"
    gid_to_mqn = {}
    for g_id, n_smi in cur.execute(qry).fetchall():
        mqn = get_MQNs(n_smi)
        if mqn:
            gid_to_mqn[g_id] = mqn
            
    # update the database with the generated MQNs
    qry = "INSERT INTO mqns VALUES " + \
          "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    for g_id in gid_to_mqn:
        qdata = (g_id, *gid_to_mqn[g_id])
        cur.execute(qry, qdata)

    print("... done")

    con.commit()
    con.close()
