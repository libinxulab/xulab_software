import sqlite3
from sqlite3 import Error

# Create a connection to the SQLite database
def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return conn

# Create a table within the SQLite database
def create_table(conn, create_table_sql):
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)

def main():
    database = "qacs.db" # Name the database

    # Name the table being created
    # Create and name columns and column types 
    sql_create_qacs_ccs_table = """ CREATE TABLE IF NOT EXISTS qacs_rt_ccs (
                                        compound_name text,
                                        exact_mz text NOT NULL,
                                        rt real,
                                        average_ccs real,
                                        ccs_calibrant text,
                                        gradient text,
                                        column_type text,
                                        notes text
                                    ); """


    # Create a database connection
    conn = create_connection(database)

    # Create tables
    if conn is not None:

        # Create projects table
        create_table(conn, sql_create_qacs_ccs_table)
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    main()