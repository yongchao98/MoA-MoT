def find_cited_lecturer():
    """
    This function provides information about the lectures Albert Einstein
    cited in his 1905 paper on Brownian motion.
    """
    # The scientist whose lectures were cited by Einstein.
    scientist_name = "Ludwig Boltzmann"
    
    # The specific work mentioned.
    work_title = "Lectures on Gas Theory"
    
    # Print the answer in a clear format.
    print("In his 1905 paper on Brownian motion, 'Über die von der molekularkinetischen Theorie der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten suspendierten Teilchen', Albert Einstein cited the work of:")
    print(f"Scientist: {scientist_name}")
    print(f"Cited work: {work_title}")

if __name__ == "__main__":
    find_cited_lecturer()