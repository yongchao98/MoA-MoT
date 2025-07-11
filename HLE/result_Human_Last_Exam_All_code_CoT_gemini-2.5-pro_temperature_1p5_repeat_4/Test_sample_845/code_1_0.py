# This script answers the user's question about the name of the value
# measured between heartbeats for HRV calculations.

# The time measured between consecutive heartbeats is the fundamental data point for HRV.
# While the user mentions PPG, the terminology is often borrowed from the gold-standard ECG (Electrocardiogram).
# In an ECG reading, the most prominent spike is the "R-wave".
# The time from one R-wave to the next is the standard measurement.
main_term = "R-R interval"

# Another common, more general term is the Inter-beat Interval.
alternative_term = "Inter-beat interval (IBI)"

print(f"The value between heartbeats used to compute HRV is called the '{main_term}' or the '{alternative_term}'.")
