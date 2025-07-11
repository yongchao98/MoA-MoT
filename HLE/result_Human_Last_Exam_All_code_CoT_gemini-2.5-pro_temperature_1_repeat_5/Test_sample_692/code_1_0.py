# The Long Count date recorded on Tres Zapotes Stela C is 7.16.6.16.18.
# This system is vigesimal (base-20) and consists of the following periods:
# 1 k'in = 1 day
# 1 winal = 20 k'in = 20 days
# 1 tun = 18 winal = 360 days
# 1 k'atun = 20 tun = 7,200 days
# 1 b'ak'tun = 20 k'atun = 144,000 days

# The date components from 7.16.6.16.18
baktun = 7
katun = 16
tun = 6
winal = 16
kin = 18

# Days per period
days_in_baktun = 144000
days_in_katun = 7200
days_in_tun = 360
days_in_winal = 20
days_in_kin = 1

# Calculate the total number of days since the Long Count calendar's start date (ca. August 11, 3114 BCE)
total_days = (baktun * days_in_baktun) + \
             (katun * days_in_katun) + \
             (tun * days_in_tun) + \
             (winal * days_in_winal) + \
             (kin * days_in_kin)

print("The earliest known recorded date from a pre-Columbian civilization is found on Tres Zapotes Stela C.")
print("The date in the Long Count calendar system is: 7.16.6.16.18")
print("\nThis corresponds to the following calculation in days:")
print(f"{baktun} * {days_in_baktun} + {katun} * {days_in_katun} + {tun} * {days_in_tun} + {winal} * {days_in_winal} + {kin} * {days_in_kin} = {total_days}")
print(f"\nThis date corresponds to September 3, 32 BCE in the Gregorian calendar.")
