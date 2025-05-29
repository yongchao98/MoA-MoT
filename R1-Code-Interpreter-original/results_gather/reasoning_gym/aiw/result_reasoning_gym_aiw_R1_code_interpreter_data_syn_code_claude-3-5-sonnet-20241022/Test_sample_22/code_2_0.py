class Circle:
    def __init__(self, name, males, females):
        self.name = name
        self.males = males
        self.females = females  # includes Matilda

def count_matilda_colleagues():
    jennifer_circle = Circle("Jennifer", 5, 3)
    joseph_circle = Circle("Joseph", 4, 4)
    
    # Count female colleagues from Jennifer's circle (excluding Matilda)
    females_from_jennifer = jennifer_circle.females - 1
    
    # Count female colleagues from Joseph's circle (excluding Matilda)
    females_from_joseph = joseph_circle.females - 1
    
    # Total female colleagues (sum because they're distinct groups)
    total_females = females_from_jennifer + females_from_joseph
    
    print(f"Female colleagues from Jennifer's circle: {females_from_jennifer}")
    print(f"Female colleagues from Joseph's circle: {females_from_joseph}")
    print(f"Total female colleagues: {total_females}")

count_matilda_colleagues()