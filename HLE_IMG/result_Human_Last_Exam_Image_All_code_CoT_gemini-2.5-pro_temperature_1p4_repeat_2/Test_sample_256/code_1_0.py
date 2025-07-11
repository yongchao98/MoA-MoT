def solve_task():
    """
    Analyzes the image to find a religious symbol not related to ancient Southwest cultures.
    """
    is_foreign_symbol_present = True
    symbol_name = "Menorah"
    description = (
        "The image is a rock art panel, likely the 'Great Gallery' in Utah, featuring Barrier Canyon Style pictographs. "
        "These are from ancient native cultures. However, one prominent figure in the upper right quadrant can be interpreted as a symbol from a different culture. "
        "This figure has a tall, ladder-like body and a head structure that strongly resembles a seven-branched candelabrum, the Jewish symbol known as a Menorah. "
        "While this is a non-standard interpretation (archaeologists consider it a native shamanic figure), the visual resemblance is notable and distinct from the other figures."
    )
    
    print(f"True or false? {is_foreign_symbol_present}")
    print(f"If true, what symbol? {symbol_name}. {description}")

solve_task()