import collections

def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    # Step 1: Define the patient's data from the case description.
    # The number 27 (patient's age) is noted here.
    patient = {
        "age": 27,
        "symptoms": ["fever", "headaches", "myalgia", "disorientation"],
        "history": "Recent camping trip to Oklahoma",
        "labs": "Elevated IgM with negative IgG Lyme serology"
    }

    # Step 2: Define the answer choices with key features.
    Disease = collections.namedtuple('Disease', ['name', 'location_match', 'notes'])
    answer_choices = {
        "A": Disease(
            name="Babesia microti",
            location_match=False,
            notes="Primarily found in the US Northeast and Upper Midwest, not Oklahoma."
        ),
        "B": Disease(
            name="Plasmodium",
            location_match=False,
            notes="Malaria is not acquired from camping in Oklahoma."
        ),
        "C": Disease(
            name="Borrelia burgdorferi",
            location_match=False,
            notes="Lyme disease is not highly endemic in Oklahoma. The Lyme test can have cross-reactivity."
        ),
        "D": Disease(
            name="Ehrlichia",
            location_match=True,
            notes="Endemic in Oklahoma (South-central US). Classic presentation with fever, headache, myalgia, and confusion. Transmitted by ticks."
        ),
        "E": Disease(
            name="Rickettsia rickettsii",
            location_match=True,
            notes="Also endemic in Oklahoma. Causes similar symptoms, but a characteristic rash (not mentioned) is common."
        )
    }

    # Step 3: Analyze and print the reasoning.
    print(f"Analyzing case for a {patient['age']}-year-old patient.")
    print(f"Key findings: {', '.join(patient['symptoms'])} after camping in Oklahoma.")
    print("-" * 30)
    print("Evaluation of possible diagnoses:")

    best_fit = None
    for key, disease in answer_choices.items():
        if disease.location_match and "Ehrlichia" in disease.name:
            best_fit = key
            print(f"-> Choice {key} ({disease.name}): Excellent Fit. {disease.notes}")
        elif disease.location_match and "Rickettsia" in disease.name:
            print(f"-> Choice {key} ({disease.name}): Good Fit. {disease.notes}")
        else:
            print(f"-> Choice {key} ({disease.name}): Poor Fit. {disease.notes}")

    print("-" * 30)
    print("Conclusion:")
    print("The combination of symptoms, recent travel to a tick-endemic area (Oklahoma),")
    print("and altered mental status is highly suggestive of Ehrlichiosis.")
    print(f"Therefore, the positive titer is expected to be for {answer_choices[best_fit].name}.")

solve_clinical_case()