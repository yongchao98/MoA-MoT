# The Long Count date from Stela 2 at Chiapa de Corzo
baktun = 7
katun = 16
tun = 3
uinal = 2
kin = 13

# Days per unit
days_per_baktun = 144000
days_per_katun = 7200
days_per_tun = 360
days_per_uinal = 20
days_per_kin = 1

# Calculate the total number of days since the Long Count epoch date
total_days = (baktun * days_per_baktun) + \
             (katun * days_per_katun) + \
             (tun * days_per_tun) + \
             (uinal * days_per_uinal) + \
             (kin * days_per_kin)

# The corresponding date in the Gregorian calendar
gregorian_date = "December 7, 36 BCE"

print("The earliest known recorded date from a pre-Columbian civilization is from Stela 2 at Chiapa de Corzo.")
print(f"Long Count Date: {baktun}.{katun}.{tun}.{uinal}.{kin}\n")
print("This represents the total number of days passed since the calendar's start date (August 11, 3114 BCE).\n")
print("Calculation:")
# The final equation with each number explicitly printed
print(f"({baktun} * {days_per_baktun}) + ({katun} * {days_per_katun}) + ({tun} * {days_per_tun}) + ({uinal} * {days_per_uinal}) + ({kin} * {days_per_kin}) = {total_days}\n")
print(f"Total days since epoch: {total_days}")
print(f"This corresponds to the Gregorian date: {gregorian_date}")