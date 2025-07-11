import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

print("Analysis of Counseling Options:")
print("=" * 35)
print("This script evaluates the best advice to give a mother whose 16-year-old son has started vaping.\n")

print("Analysis of Statement I:")
print("   'Vaping is a better option for her son than cigarettes. It is ok for him to continue vaping, as long as he is trying to cut down.'")
print("   - Verdict: Incorrect. While vaping is less harmful than smoking, it is not harmless. For an adolescent, the primary goal is cessation from all nicotine products due to risks to brain development and the potential for long-term addiction. Normalizing vaping is not recommended.\n")

print("Analysis of Statement II:")
print("   'It would be good for her son to start using nicotine patches, gum, or lozenges instead of vaping.'")
print("   - Verdict: Correct. Nicotine Replacement Therapy (NRT) is a medically approved, first-line strategy to help people quit nicotine. It allows the user to manage withdrawal symptoms safely without the harmful inhalants from vapes or cigarettes.\n")

print("Analysis of Statement III:")
print("   'Vaping’s risks and benefits are well understood in adults, but not in children, so her son should not vape at all.'")
print("   - Verdict: Correct. This is a critical point. The long-term health effects of vaping are unknown, especially on the developing adolescent body and brain. This statement correctly frames why her son's situation is different from her own (an adult smoker using vapes for harm reduction) and supports the goal of complete cessation.\n")

print("Analysis of Statement IV:")
print("   'Vaping’s risks and benefits remain poorly understood in children, but it has been shown vaping has clear benefits over cigarettes in children.'")
print("   - Verdict: Incorrect. This phrasing is misleading. Labeling anything that is still harmful and addictive as having 'clear benefits' for a child is dangerous and inappropriate. The comparison should not be used to justify adolescent vaping.\n")

print("Analysis of Statement V:")
print("   'Consider initiating bupropion and varenicline depending on her son’s needs.'")
print("   - Verdict: Plausible, but a secondary consideration. These are prescription medications for smoking cessation. While they can be used off-label in adolescents, they are typically considered second-line options after behavioral counseling and NRT (Statement II) have been attempted. They are not usually the starting point for counseling.\n")

print("=" * 35)
print("Conclusion:")
print("The best counseling strategy combines explaining the rationale for quitting (the unknown risks for adolescents) and providing a concrete, safe, first-line treatment option. Therefore, the combination of statements II and III is the most appropriate.")

print("\nThe recommended options for the final answer are:")
print("II")
print("III")

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)

# The letter for the combination of II and III is J.
final_answer = "J"
print(f'<<<J>>>')