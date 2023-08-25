import sqlite3
import pandas as pd

# Connect to the initialized database where data will be entered
cxn = sqlite3.connect("qacs.db") 
wb = pd.read_excel("QAC Database_v3.xlsx", sheet_name = "qacs_rt_ccs") # Insert path to Excel spreadsheet containing data
wb.to_sql(name="qacs_rt_ccs", con=cxn, if_exists="replace", index=True)

cxn.commit()
cxn.close()