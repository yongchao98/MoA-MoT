import decimal

decimal.getcontext().prec = 50

C = decimal.Decimal(25200)
M = decimal.Decimal(3628800)

log10_M = M.log10()
log10_N = C.log10() + M * log10_M

integer_part = int(log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR))
fractional_part = log10_N - integer_part

leading_digits_factor = decimal.Decimal(10)**(fractional_part)

final_string = '{:.49f}'.format(leading_digits_factor)
first_40_digits = final_string.replace('.', '')[:40]

print(first_40_digits)