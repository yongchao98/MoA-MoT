# This script will define and print the name of the value computed between heartbeats for HRV analysis.

# The value is the time interval between consecutive heartbeats.
# In electrocardiography (ECG), this is measured from the peak of one QRS complex (the R wave) to the next.
# This term is widely used even when data is collected from other sources like PPG.
value_name = "R-R interval"

print(f"The value computed between heartbeats to measure HRV is called the: {value_name}")
