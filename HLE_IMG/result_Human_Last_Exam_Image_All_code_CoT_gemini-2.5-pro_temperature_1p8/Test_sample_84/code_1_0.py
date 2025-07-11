# The plan is to identify the city based on the distinct heraldic elements shown in the image.
# Step 1: Analyze the shield. It is divided horizontally (per fess).
# The top field is blue (azure) and contains a yellow crescent moon (a crescent Or).
# The bottom field is yellow (or) and contains a black crab (a crab Sable).
# Step 2: Search for a municipal coat of arms with this combination.
# The combination of a crab and a crescent moon on a blue and yellow field is very specific.
# Step 3: The search points to the coat of arms of Ochakiv, a city in southern Ukraine.
# Its coat of arms perfectly matches the description of the shield in the image.
# Step 4: Conclude that the city in question is Ochakiv.
# The following code will print the name of the identified city.

def find_city_by_coat_of_arms():
    """
    Identifies and prints the name of the city with a coat of arms
    matching the visual description.
    """
    # Key elements from the shield
    top_color = "Blue"
    top_charge = "Crescent Moon"
    bottom_color = "Gold"
    bottom_charge = "Crab"

    # The city identified based on these heraldic elements
    city_name = "Ochakiv"

    print(f"Based on the elements of the shield ({top_color} field with a {top_charge}, {bottom_color} field with a {bottom_charge}), the city is {city_name}.")

find_city_by_coat_of_arms()