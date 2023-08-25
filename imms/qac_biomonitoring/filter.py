import pandas as pd

def filter_data(output_file):
    df = pd.read_excel(output_file) # Read the input processed file into a pandas dataframe

    # Filter rows where sample_type = sample
    sample_rows = df[df["sample_type"] == "sample"]

    # Create an empty list to store unique features
    unique_features = []

    # Iterate over sample rows
    for _, sample_row in sample_rows.iterrows():
        sample_mz = sample_row["mz"]
        sample_rt = sample_row["rt"]
        sample_peak_area = sample_row["peak_area"]
        
        # Filter control rows with matching m/z and rt within +/- 0.1
        control_rows = df[(df["sample_type"] == "control") &
                          (df["mz"] == sample_mz) &
                          (abs(df["rt"] - sample_rt) <= 0.1)] # Adjust rt tolerance 
        
        # Check if there are any corresponding control peaks
        if len(control_rows) == 0:
            unique_features.append(sample_row)  # No corresponding control peak, add as unique feature
        else:
            control_peak_area = control_rows["peak_area"].max()
            if sample_peak_area > 1.5 * control_peak_area: # Check if sample peak area is 1.5 times greater than control peak area
                unique_features.append(sample_row)
               
    # Create a new dataframe from the unique features list
    output_df = pd.DataFrame(unique_features)

    # Write the output DataFrame to an Excel file
    filtered_file = output_file[:-15] + "_filtered.xlsx"
    output_df.to_excel(filtered_file, index=False)

    return filtered_file