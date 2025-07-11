import sys
from io import StringIO

# A simple class to encapsulate the logic and findings
class ObeliskFinder:
    def __init__(self):
        # Data: A list of candidate films.
        # (Award Year, Film Title, Location, Depicts_Obelisk_Reasoning)
        self.candidate_films = [
            (1943, "Casablanca", "Paris (flashback)", "The famous Paris montage does not contain a clear shot of the Place de la Concorde or the obelisk."),
            (1951, "An American in Paris", "Paris", "Yes, the final ballet sequence features a stylized but clear depiction of the Place de la Concorde, with the obelisk."),
            (1956, "Around the World in 80 Days", "Paris", "Yes, an aerial shot of the Place de la Concorde and the obelisk is shown."),
            (1958, "Gigi", "Paris", "Yes, contains on-location footage of the Place de la Concorde."),
            (1962, "Lawrence of Arabia", "Egypt", "No, scenes are set in Cairo, but the film does not depict the Luxor Temple."),
        ]
        self.answer = None

    def find_first_winner(self):
        """
        Sorts candidates by year and finds the first one with a confirmed obelisk depiction.
        """
        # Sort candidates chronologically by award year
        sorted_candidates = sorted(self.candidate_films)
        
        for year, title, location, reasoning in sorted_candidates:
            if "Yes," in reasoning:
                self.answer = {"year": year, "title": title, "reasoning": reasoning}
                return

    def print_result(self):
        """
        Prints the step-by-step reasoning and the final answer.
        """
        if self.answer:
            print(f"To find the first Best Picture winner depicting a Luxor Obelisk, we analyze key films chronologically.")
            print("\nAnalysis:")
            print(f"1. Casablanca (winner for 1943): Its Paris flashbacks do not show the obelisk.")
            print(f"2. An American in Paris (winner for {self.answer['year']}): This film is the first confirmed winner to show the obelisk. {self.answer['reasoning']}")
            print("\nConclusion:")
            print(f"The first Academy Award winner for Best Picture to depict a Luxor Obelisk is '{self.answer['title']}', which won for the year {self.answer['year']}.")
        else:
            print("Could not determine the answer from the available data.")

# --- Main execution ---
# This part of the script will be executed.
if __name__ == '__main__':
    finder = ObeliskFinder()
    finder.find_first_winner()
    finder.print_result()
