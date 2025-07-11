def explain_process():
  """
  This function explains the biological process shown in the provided image sequence.
  """
  title = "Process Depicted: Hyphal Branching"
  explanation_lines = [
      "The image sequence (L, M, N) illustrates the process of hyphal branching in a filamentous fungus.",
      "",
      "Here is a step-by-step breakdown:",
      "1. Panel L: Shows a mature, septate (containing cross-walls) fungal hypha. The structure is relatively straight and uniform. The arrowhead points to a septum.",
      "2. Panel M: A small bulge begins to form on the side of the hypha, near the septum. This is the initiation of a new branch, as highlighted by the arrowhead.",
      "3. Panel N: The bulge has grown outwards to form a new, distinct lateral branch. This process allows the fungus to expand its network (mycelium) to search for nutrients."
  ]
  
  print(title)
  print("-" * len(title))
  for line in explanation_lines:
    print(line)

explain_process()