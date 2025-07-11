def find_cited_lectures():
    """
    Provides information about citations in Einstein's 1905 papers.
    """
    
    paper_in_question = "1905 paper on Brownian motion"
    cited_scientist_in_brownian_paper = "None"
    
    other_famous_paper = "1905 paper on the photoelectric effect"
    cited_scientist_in_other_paper = "Ludwig Boltzmann"
    cited_work = "Vorlesungen über Gastheorie (Lectures on Gas Theory)"
    
    print(f"Strictly speaking, in his {paper_in_question}, Albert Einstein cited {cited_scientist_in_brownian_paper}.")
    print("The paper is famous for not referencing prior work.")
    print("\nHowever, in his {other_famous_paper} from the same year, Einstein did cite lectures from a key figure:")
    print(f"\nScientist: {cited_scientist_in_other_paper}")
    print(f"Work Cited: {cited_work}")
    print("\nBoltzmann's work was fundamental to the ideas in the Brownian motion paper, even though he was not formally cited within it.")

find_cited_lectures()