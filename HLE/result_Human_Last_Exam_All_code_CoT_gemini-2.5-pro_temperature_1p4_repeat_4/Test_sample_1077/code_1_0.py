def find_philosopher():
    """
    Identifies the philosopher who wrote the specified letter to Heidegger.
    """
    philosopher_name = "Karl Löwith"
    
    explanation = (
        "The letter was written by the German philosopher Karl Löwith in 1936.\n"
        "Löwith, a former student of Heidegger, wrote this letter after a "
        "disappointing meeting with him in Rome.\n"
        "He expressed his deep disillusionment with Heidegger's character and his "
        "unrepentant support for Nazism, which stood in stark contrast to the "
        "profound philosophical insights Löwith admired in 'Being and Time'."
    )
    
    print(f"The philosopher who wrote the letter to Heidegger was: {philosopher_name}")
    print("\n--- Context ---")
    print(explanation)

if __name__ == "__main__":
    find_philosopher()