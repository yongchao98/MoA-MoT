import sys

def solve():
  """
  This function solves the nyala re-identification task.
  
  To identify the correct nyala, we look for unique, unchanging features.
  The most reliable features are the horn shape and the pattern of white spots on the flank.

  1.  **Target Nyala's Features:**
      -   **Horns:** Long, lyre-shaped, with a specific curve on the right horn.
      -   **Spot Pattern:** A distinct cluster of 3-4 white spots is visible on the flank, behind the shoulder. This pattern is unique to the individual.

  2.  **Comparison with Options:**
      -   **A:** Horns are smaller; likely a younger, different animal.
      -   **B:** Features are similar, a good candidate.
      -   **C:** Image quality is poor, but horn shape seems to match.
      -   **D:** The horn shape and, critically, the unique pattern of white spots on the flank are an identical match to the target.
      -   **E:** Has a white band/collar on its leg that the target lacks. Definitely a different animal.
      -   **F:** Partially obscured, preventing a conclusive comparison.

  3.  **Final Decision:**
      Image D provides the clearest and most definitive match to the target animal, based on the identical horn shape and the unique "fingerprint" of its spot pattern.
  """
  # The correct image is D.
  answer = 'D'
  print(f"{answer}")

solve()