import webbrowser

def manual_counter():
    """
    A tool to help you manually count specific events in a video.
    """
    video_url = "https://www.youtube.com/watch?v=F553ma0hoxc"
    
    print("--- Pas de Balanchine Counter for Svetlana Gileva's Solo ---")
    print("\nThis script will help you count the number of 'pas de Balanchines'.")
    print("\nStep 1: Get Ready")
    print("A 'pas de Balanchine' is a specific ballet step. It looks like a 'pas de bourrée piqué' followed by an arabesque that lowers into a 'plié'.")
    print("It has a distinct 'step-step-step-lunge' rhythm.")
    
    print("\nStep 2: Watch the Video")
    print(f"The solo scene can be found here: {video_url}")
    print("Your web browser will now open to the video.")
    
    try:
        webbrowser.open(video_url)
    except Exception:
        print("Could not open web browser automatically. Please copy the link above.")

    print("\nStep 3: Count the Steps")
    print("Position the video window next to this terminal.")
    print("Press the [Enter] key each time Svetlana Gileva performs a pas de Balanchine.")
    print("Type 'q' and then press [Enter] when you are finished or want to quit.")
    print("-" * 60)

    count = 0
    while True:
        user_input = input(f"Current count: {count}. Press Enter to count the next one, or type 'q' to finish: ")
        if user_input.lower() == 'q':
            break
        count += 1
        print(f"Step number {count} counted.")

    print("\n==============================================")
    print("             FINAL RESULT")
    print("==============================================")
    # The "equation" is a simple summation of 1s.
    # For a final count of N, the equation is 1 + 1 + ... (N times).
    if count > 0:
        equation_str = " + ".join(["1"] * count)
        print(f"Equation: {equation_str} = {count}")
    else:
        print("Equation: 0")
    
    print(f"\nSvetlana Gileva performed a total of {count} pas de Balanchines.")

if __name__ == "__main__":
    manual_counter()
