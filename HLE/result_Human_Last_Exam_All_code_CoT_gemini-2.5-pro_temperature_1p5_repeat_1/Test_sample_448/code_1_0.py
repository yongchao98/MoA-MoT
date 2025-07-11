# The rarest stable noble gas in Earth's atmosphere is Xenon.
gas_name = "Xenon"
# Its abundance in the atmosphere is approximately 0.087 parts per million (ppm).
abundance_ppm = 0.087

# Convert parts per million (ppm) to percentage.
# 1 ppm = 0.0001%, so we divide the ppm value by 10,000.
denominator = 10000
percentage = abundance_ppm / denominator

# Print the final result in the format of an equation.
print(f"The rarest stable noble gas is {gas_name}.")
print("Its percentage in Earth's atmosphere is calculated as follows:")
print(f"{abundance_ppm} ppm / {denominator} = {percentage}%")