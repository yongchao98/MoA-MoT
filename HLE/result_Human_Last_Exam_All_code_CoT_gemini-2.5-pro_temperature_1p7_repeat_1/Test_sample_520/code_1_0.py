import string

def solve():
    """
    Solves the trivia question about Andre Vernault by constructing an
    "equation" from the letters of his profession, as requested.
    """
    # The answer to the historical question, based on my knowledge base.
    profession = "smuggler"

    # Create a list to hold the numeric value of each letter (a=1, b=2, etc.).
    numbers = []
    for char in profession:
        # 'a' is 1, 'b' is 2, and so on.
        number = string.ascii_lowercase.index(char) + 1
        numbers.append(str(number))

    # Format the numbers and the final word into an "equation".
    equation_str = " + ".join(numbers)

    # Print the final equation as requested.
    print("Here is the equation to reveal the profession:")
    print(f"{equation_str} = {profession}")
    print("-" * 30)
    print(f"Based on historical records, Andre Vernault's concealed true profession was: {profession}")

solve()