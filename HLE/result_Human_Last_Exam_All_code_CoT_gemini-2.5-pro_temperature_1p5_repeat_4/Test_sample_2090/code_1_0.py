# The IAU constellation boundaries are defined for the B1875.0 epoch.
# The boundary where the Vernal Equinox precessed into Pisces is a line
# of constant Right Ascension. These are the two vertices defining the
# segment that crosses the celestial equator (Declination = 0).

# Coordinates for the first point (with the lower declination)
ra_h = 23
ra_m = 58
ra_s = 22
dec_d1 = -2
dec_m1 = 49
dec_s1 = 0

# Coordinates for the second point (with the higher declination)
# RA is the same for both points as they form a line of constant RA.
dec_d2 = 1
dec_m2 = 8
dec_s2 = 0

# --- Calculations ---

# Convert Declinations to decimal degrees.
# The formula is: Decimal Degrees = Degrees + (Minutes / 60) + (Seconds / 3600)
dec_decimal1 = dec_d1 - (dec_m1 / 60) - (dec_s1 / 3600)
dec_decimal2 = dec_d2 + (dec_m2 / 60) + (dec_s2 / 3600)

# --- Formatting ---

# Format the Right Ascension string "XX YY ZZ"
ra_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}"

# Format the Declination strings "AA.BB" with a leading zero for the integer part.
# The format needs to handle negative signs correctly.
def format_decimal_declination(dec_val):
    """Formats declination to a signed string like '-02.82' or '01.13'."""
    sign = "-" if dec_val < 0 else ""
    # We round the absolute value to two decimal places to handle precision.
    abs_rounded_val = round(abs(dec_val), 2)
    integer_part = int(abs_rounded_val)
    # The decimal part is the remainder after getting the integer.
    decimal_part = int(abs_rounded_val * 100) % 100
    # The problem implies we do not want a '+' sign for positive values.
    return f"{sign}{integer_part:02d}.{decimal_part:02d}"

dec_str1 = format_decimal_declination(dec_decimal1)
dec_str2 = format_decimal_declination(dec_decimal2)

# Assemble the final string, with the lower point first.
final_output = f"{ra_str}, {dec_str1}; {ra_str}, {dec_str2}"

# Print the final result
print(final_output)
