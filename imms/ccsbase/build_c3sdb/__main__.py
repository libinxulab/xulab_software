"""
    build_c3sdb/__main__.py
    Dylan H. Ross

    this script coordinates the high-level steps in the database building process 
"""

import os

from .initialize_db_tables import main as initialize_db_tables_main
from .fill_db_from_src import main as fill_db_from_src_main
#from .fetch_smiles import main as fetch_smiles_main
#from .generate_mqns import main as generate_mqns_main
#from .label_chem_class import main as label_chem_class_main


def main():
    """ main database build sequence """

    # remove C3S.db if it exists
    if os.path.isfile('C3S.db'):
        os.remove('C3S.db')

    # initialize the database with empty data tables
    initialize_db_tables_main()

    # add data from cleaned datasets into the database
    fill_db_from_src_main()

    # try to fill in missing SMILES structures
    #fetch_smiles_main()

    # add in MQNs for all available SMILES structures
    #generate_mqns_main() 

    # label each compound with an approximate chemical classification
    #label_chem_class_main()


if __name__ == '__main__':
    main()
