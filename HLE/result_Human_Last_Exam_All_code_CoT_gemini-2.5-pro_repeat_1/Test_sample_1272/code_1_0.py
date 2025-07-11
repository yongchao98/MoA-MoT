def find_defeat_location():
    """
    Identifies and prints the location of Ming General Mu Sheng's first major defeat.
    """
    general = "Ming General Mu Sheng"
    event = "the Battle of Tốt Động–Chúc Động"
    year = 1426
    
    # While Wang Tong was the field commander, this was the first catastrophic defeat
    # for the Ming army during the Lam Sơn uprising, which Mu Sheng was part of.
    ancient_location = "Tốt Động and Chúc Động"
    modern_district = "Chương Mỹ District"
    modern_province = "Hanoi"

    print(f"The first major defeat experienced by {general}'s forces during the Lam Sơn uprising occurred in {year} at {event}.")
    print(f"The battle took place at {ancient_location}, which is now located in {modern_district}.")
    print(f"This area is part of the modern-day province-level municipality of: {modern_province}")

find_defeat_location()