def find_cited_lecturer():
    # This script identifies the scientist whose lectures Einstein's 1905 Brownian motion paper relied on.
    # The work was not a direct citation of lectures Einstein attended, but rather of published lectures he studied.
    
    scientist_name = "Ludwig Boltzmann"
    work_title = "Lectures on Gas Theory"
    year = 1905
    paper_topic = "Brownian motion"

    # Construct and print the answer
    print(f"In his {year} paper on {paper_topic}, Albert Einstein's work was fundamentally based on the principles from the published lectures of {scientist_name}, specifically his work titled '{work_title}'.")

if __name__ == "__main__":
    find_cited_lecturer()