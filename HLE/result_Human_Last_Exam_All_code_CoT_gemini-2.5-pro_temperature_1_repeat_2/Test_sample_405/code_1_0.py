import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a new in-memory text stream
captured_output = io.StringIO()
# Redirect stdout to the new stream
sys.stdout = captured_output

# Phosphorus has several allotropes, each with a distinct color.
# We will count the number of unique colors from the most common and stable forms.

# 1. White Phosphorus: This allotrope is a waxy solid that is white,
#    but can appear yellowish due to impurities or conversion to red phosphorus.
#    We will count this as 1 color: White.
color_1 = "White"
num_1 = 1

# 2. Red Phosphorus: This is an amorphous form whose color can vary
#    from orange to a deep red. We will count this as 1 color: Red.
color_2 = "Red"
num_2 = 1

# 3. Violet Phosphorus: Also known as Hittorf's metallic phosphorus,
#    this is a well-defined crystalline form that is deep violet.
#    We will count this as 1 color: Violet.
color_3 = "Violet"
num_3 = 1

# 4. Black Phosphorus: This is the most stable form, with a layered
#    structure similar to graphite. It is black and has a metallic luster.
#    We will count this as 1 color: Black.
color_4 = "Black"
num_4 = 1

# Now, we sum the counts for each distinct color.
total_colors = num_1 + num_2 + num_3 + num_4

print("The distinct colors of phosphorus allotropes are:")
print(f"- {color_1}")
print(f"- {color_2}")
print(f"- {color_3}")
print(f"- {color_4}")
print("\nTo find the total number of colors, we can represent it as an equation:")
# The instruction requires outputting each number in the final equation.
print(f"{num_1} ({color_1}) + {num_2} ({color_2}) + {num_3} ({color_3}) + {num_4} ({color_4}) = {total_colors}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output_str = captured_output.getvalue()

# Print the final result to the user
print(output_str)

# The final answer in the required format
print(f"<<<{total_colors}>>>")