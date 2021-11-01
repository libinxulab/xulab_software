"""
    initialize_db_tables.py
    Dylan H. Ross

    Initializes a new database with empty tables
"""

from sqlite3 import connect


# define the table schemas: master, mqns, predicted

# main table for combining the measurement data
master_schema = """
CREATE TABLE master (
    -- global unique string identifier
    g_id TEXT UNIQUE NOT NULL,
    -- compound name
    name TEXT NOT NULL,
    -- MS adduct
    adduct TEXT NOT NULL,
    -- mass, z, m/z and CCS
    mass REAL NOT NULL,
    z INTEGER NOT NULL,
    mz REAL NOT NULL,
    ccs REAL NOT NULL,
    -- neutral smiles structure
    smi TEXT,
    -- (rough) chemical class label
    chem_class_label TEXT,
    -- tag referencing which dataset the value is from
    src_tag TEXT NOT NULL,
    -- CCS type (DT, TW, ...)
    ccs_type TEXT NOT NULL,
    -- describe method used for CCS measurement (e.g. stepped-field, calibrated with polyalanine)
    ccs_method TEXT NOT NULL
);
"""

# table for molecular descriptors (MQNs)
mqns_schema = """
CREATE TABLE mqns (
    -- global unique integer identifier (same as in master)
    g_id TEXT UNIQUE NOT NULL,
    -- atom counts (12)
    c INTEGER NOT NULL, f INTEGER NOT NULL, cl INTEGER NOT NULL, br INTEGER NOT NULL,
    i INTEGER NOT NULL, s INTEGER NOT NULL, p INTEGER NOT NULL, an INTEGER NOT NULL,
    cn INTEGER NOT NULL, ao INTEGER NOT NULL, co INTEGER NOT NULL, hac INTEGER NOT NULL,
    -- polarity counts (6)
    hbam INTEGER NOT NULL, hba INTEGER NOT NULL, hbdm INTEGER NOT NULL,
    hbd INTEGER NOT NULL, neg INTEGER NOT NULL, pos INTEGER NOT NULL,
    -- bond counts (7)
    asb INTEGER NOT NULL, adb INTEGER NOT NULL, atb INTEGER NOT NULL, csb INTEGER NOT NULL,
    cdb INTEGER NOT NULL, ctb INTEGER NOT NULL, rbc INTEGER NOT NULL,
    -- topology counts (17)
    asv INTEGER NOT NULL, adv INTEGER NOT NULL, atv INTEGER NOT NULL, aqv INTEGER NOT NULL,
    cdv INTEGER NOT NULL, ctv INTEGER NOT NULL, cqv INTEGER NOT NULL, r3 INTEGER NOT NULL,
    r4 INTEGER NOT NULL, r5 INTEGER NOT NULL, r6 INTEGER NOT NULL, r7 INTEGER NOT NULL,
    r8 INTEGER NOT NULL, r9 INTEGER NOT NULL, rg10 INTEGER NOT NULL, afr INTEGER NOT NULL,
    bfr INTEGER NOT NULL
);
"""

# table for predicted CCS values
predicted_schema = """
CREATE TABLE predicted (
    -- global unique integer identifier (same as in master, increment past highest value in master for new predictions)
    g_id TEXT UNIQUE NOT NULL,
    -- compound name
    name TEXT NOT NULL,
    -- MS adduct
    adduct TEXT NOT NULL,
    -- m/z and CCS (predicted)
    mz REAL NOT NULL,
    pred_ccs REAL NOT NULL,
    -- neutral smiles structure
    smi TEXT NOT NULL,
    -- class label from untargeted classification step (if used)
    class_label INTEGER,
    -- CCS error relative to reference value in master table (if available)
    pred_error REAL,
    -- timestamp (YYMMDDHHmmss) of when prediction was added to the database, NULL for original reference data
    t_stamp INTEGER
)
;"""


def create_tables(cursor):
    """ creates all of the database tables """
    cursor.execute(master_schema)
    cursor.execute(mqns_schema)
    cursor.execute(predicted_schema)


def main():
    """ main execution sequence """

    # initialize the database connection
    con = connect('C3S.db')
    cur = con.cursor()

    # create the tables
    create_tables(cur)

    # commit changes and close DB connection
    con.commit()
    con.close()


