import sys

def evaluate_mausoleum_claims():
    """
    Evaluates historical claims about the Mausoleum of Augustus to find the one
    with the most backing from archaeological historians.
    """

    # Each option is broken into individual claims, which are scored for accuracy.
    # Score key: 1.0 = Correct, 0.5 = Partially correct/Debatable, 0.0 = Incorrect/Subjective
    options = {
        'A': {
            'claims': {
                "Is a demonstration of the emperor's great power": 1.0,
                "Power peaked in his fifth decade": 0.5,
                "Power waned as he grew older": 0.5
            },
            'rationale': "A plausible interpretation of political history, but the specific timeline of his power's peak and wane is debatable and not a primary archaeological conclusion."
        },
        'B': {
            'claims': {
                "Could be likened to the tomb of Mausolus": 1.0,
                "Was smaller than the tombs of the kings of Namibia": 0.0
            },
            'rationale': "The comparison to Mausolus is correct (it's the origin of the word 'mausoleum'), but the comparison to Namibian tombs is unsubstantiated and likely false."
        },
        'C': {
            'claims': {
                "Strabo claimed it had evergreen trees on top": 1.0,
                "Half amphoras were placed on the sides for libations": 0.0,
                "Mounds were graves for his relatives and friends": 1.0,
                "An elevated place where Livia's corpse was burnt was nearby": 1.0
            },
            'rationale': "Mostly correct based on the primary source Strabo and archaeology (the ustrinum for cremation), but Strabo does not mention the half amphoras, a specific detail that makes the statement inaccurate."
        },
        'D': {
            'claims': {
                "Tacitus claimed it cost millions in silver denarii": 0.0,
                "Was designed by slaves from Judea": 0.0
            },
            'rationale': "This statement is anachronistic and unsubstantiated. The Mausoleum was built before the major enslavement of Judeans following the Jewish-Roman War of 70 AD."
        },
        'E': {
            'claims': {
                "Official designation was Tumulus Iulorium": 0.5,
                "Some called it Stultitia Augusti ('Augustus's Folly')": 0.0
            },
            'rationale': "While a plausible Latin name, 'Mausoleum Augusti' was its known name. There is no historical evidence for the nickname 'Stultitia Augusti'."
        },
        'F': {
            'claims': {
                "Adopted an architectural style for self-glorification of eastern dynasts": 1.0,
                "The full text of the Res Gestae Divi Augusti was inscribed on it": 1.0
            },
            'rationale': "Both claims are core, well-established facts. The link to Hellenistic models is a foundational interpretation, and the presence of the Res Gestae on bronze pillars at the entrance is confirmed by the text of the Res Gestae itself."
        },
        'G': {
            'claims': {
                "Was taller than the Mausoleum of Mausolos": 0.5
            },
            'rationale': "This is a debatable claim. Reconstructions of both vary, and it is not universally agreed that Augustus's Mausoleum was definitively taller. It was of a comparable, massive scale."
        },
        'H': {
            'claims': {
                "The mausoleum of Hadrian has a clearer architectural language": 0.0
            },
            'rationale': "This is a subjective statement of architectural criticism, not a verifiable historical or archaeological fact."
        },
        'I': {
            'claims': {
                "Was capped by a conical roof and a bronze statue of Augustus": 1.0,
                "Included a chapel built to the Archangel Michael": 0.0
            },
            'rationale': "This statement incorrectly mixes the original Roman structure with a much later medieval addition. The chapel was built in the 12th century when the ruin was converted into a fortress."
        }
    }

    best_option = None
    max_score = -1.0

    print("Analyzing each option by scoring its claims...\n")

    for option_key, data in sorted(options.items()):
        claims = data['claims']
        # The line below suppresses a false positive from mypy when running on Python 3.12+
        # It's not essential for the logic but is good practice.
        if sys.version_info >= (3, 12) and not claims: # type: ignore
            continue
        scores = list(claims.values())
        avg_score = sum(scores) / len(scores)

        # Outputting the 'equation' for each option's score
        # The format string below correctly formats numbers like 1.0 as '1'.
        score_numbers = [f"{score:g}" for score in scores]
        equation = f"({ ' + '.join(score_numbers) }) / {len(scores)}"

        print(f"Option {option_key}:")
        print(f"  Rationale: {data['rationale']}")
        print(f"  Score Calculation: {equation} = {avg_score:.2f}\n")

        if avg_score > max_score:
            max_score = avg_score
            best_option = option_key

    print("---" * 10)
    print(f"Conclusion: Option '{best_option}' has the highest score ({max_score:.2f}), indicating it has the most backing from archaeological and historical sources.")
    print(f"<<<{best_option}>>>")


if __name__ == '__main__':
    evaluate_mausoleum_claims()