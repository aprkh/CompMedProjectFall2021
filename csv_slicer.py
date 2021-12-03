import pandas as pd

# this little helper file is just for splicing out columns of interest
# since running my code repeatedly with the whole giant CSVs takes forever

df = pd.read_csv('ctrl_vs_case.csv')

# isolate the columns I want
df = df[['Participant_ID', 'CtrlVsCase_Classifier', 'DEAF1', 'CCR1', 'STAT4', 'BEAN1', 'TACR3', 'GPAA1', 'NDUFS7', 'CC2D1A', 'GRM1', 'GAMT']]

# output the csv
df.to_csv('case_10.csv')
