import re

def check_correctness_of_biology_riddle(llm_answer_text: str) -> str:
    """
    Checks the correctness of the answer to the biology riddle about protein secretion.

    The function deciphers the riddle's clues to determine the correct start and end
    locations and compares this to the provided answer.
    """

    # 1. Define the biological knowledge derived from the riddle's clues.
    # "ribonucleoprotein particle" (SRP) meets "nascent chain" (on a ribosome) -> cytosol
    meeting_place = "cytosol"
    
    # "sugar" (glycosylation), "rough" (RER), "on my way" (secretory pathway)
    # -> final destination is outside the cell.
    destination = "extracellular space"

    # 2. Define the multiple-choice options.
    options = {
        "A": ("cytosol", "extracellular space"),
        "B": ("ribosome", "proteasome"),
        "C": ("membrane", "nucleus"),
        "D": ("Golgi", "mitochondrion")
    }

    # 3. Determine the correct option key based on the analysis.
    correct_option_key = None
    for key, (start, end) in options.items():
        if start == meeting_place and end == destination:
            correct_option_key = key
            break
    
    if correct_option_key is None:
        # This should not be reached if the question and options are valid.
        return "Error in checking logic: Could not find a matching option for the deduced path."

    # 4. Extract the final answer from the LLM's response.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "The provided answer does not contain a final choice in the format <<<X>>>."
    
    llm_answer_key = match.group(1)

    # 5. Compare the LLM's answer with the correct answer and return the result.
    if llm_answer_key == correct_option_key:
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect.\n"
            f"1. **Meeting Place Analysis**: The riddle describes a 'ribonucleoprotein particle' (Signal Recognition Particle, SRP) meeting a 'nascent chain' (a new protein on a ribosome). This interaction occurs in the **cytosol**.\n"
            f"2. **Destination Analysis**: The clues 'you really need some sugar' (glycosylation), 'it seems somewhat rough' (Rough Endoplasmic Reticulum), and 'I need to be on my way' (secretory pathway) all point to a protein being secreted from the cell. The final destination is the **extracellular space**.\n"
            f"3. **Conclusion**: The correct path is from the 'cytosol' to the 'extracellular space', which corresponds to option {correct_option_key}.\n"
            f"The provided answer was {llm_answer_key}, which represents the path '{' to '.join(options[llm_answer_key])}'. This is incorrect."
        )
        return reason

# The final answer provided by the LLM is passed to the checking function.
llm_response = """
Analyse der Frage und der Kandidatenantworten:

1.  **Analyse der Rätselfrage**: Die Frage beschreibt einen zellbiologischen Prozess in Form eines Dialogs.
    *   **"Oh, ich kenne dich", sagt das Ribonukleoproteinpartikel zur naszierenden Kette, als sie sich treffen.** Dies bezieht sich auf das Signalerkennungspartikel (SRP), ein Ribonukleoprotein, das eine neu synthetisierte Proteinkette ("naszierende Kette") an einem Ribosom erkennt. Dieses Treffen findet im **Zytosol** statt, wo die Proteinsynthese beginnt.
    *   **"Halte dort für eine Minute an. Lass mich dich hineinführen; du brauchst wirklich etwas Zucker."** Das SRP bewirkt eine vorübergehende Unterbrechung der Proteinsynthese ("Pause"). Es führt dann den gesamten Komplex zum endoplasmatischen Retikulum (ER). Der Hinweis auf "Zucker" bezieht sich auf die Glykosylierung, eine Modifikation, die im ER beginnt.
    *   **"Es scheint etwas rau zu sein."** Dies ist ein direkter Hinweis auf das **raue endoplasmatische Retikulum** (RER), das aufgrund der Ribosomen auf seiner Oberfläche "rau" erscheint.
    *   **"Ich schätze, das ist ein Abschied; ich muss weiter", antwortet die Kette.** Nachdem das SRP den Komplex zum RER gebracht hat, wird es freigesetzt ("Abschied"). Die Proteinkette tritt in den Sekretionsweg ein ("ich muss weiter"). Ein primäres Endziel für Proteine auf diesem Weg ist die Sekretion aus der Zelle.

2.  **Bestimmung des Weges**:
    *   **Treffpunkt**: Das SRP trifft die naszierende Kette im **Zytosol**.
    *   **Zielort**: Die Kette wird in das RER transportiert, durch den Golgi-Apparat weiterverarbeitet und schließlich aus der Zelle in den **extrazellulären Raum** ausgeschieden.

3.  **Bewertung der Optionen**:
    *   **A) Zytosol zum extrazellulären Raum**: Dies beschreibt den Startpunkt und den endgültigen Zielort des im Rätsel beschriebenen sekretorischen Weges korrekt.
    *   **B) Ribosom zum Proteasom**: Das Proteasom ist für den Proteinabbau zuständig, nicht für die Sekretion.
    *   **C) Membran zum Zellkern**: Das Treffen findet im Zytosol statt, nicht an einer Membran, und der Weg führt nicht zum Zellkern.
    *   **D) Golgi zum Mitochondrium**: Das Treffen findet vor dem Golgi statt, und dies ist nicht der Weg zum Mitochondrium.

4.  **Analyse der Kandidatenantworten**: Fast alle Kandidaten, die eine Begründung liefern (Antworten 1-12 und 15), beschreiben den korrekten biologischen Prozess und identifizieren den Weg vom Zytosol in den extrazellulären Raum. Viele dieser Antworten geben jedoch am Ende einen falschen Buchstaben an (z. B. D, C oder B), was auf einen systematischen Fehler bei der Zuordnung der korrekten Begründung zur richtigen Antwortoption hindeutet. Die Begründungen selbst unterstützen jedoch durchweg die Option A. Die Antworten 5, 8, 9, 11 und 14 kommen zum korrekten Schluss.

**Schlussfolgerung**: Basierend auf der schrittweisen Analyse des Rätsels und der überwältigenden Übereinstimmung in den Begründungen der Kandidaten ist der korrekte Weg vom Zytosol in den extrazellulären Raum.

<<<A>>>
"""

result = check_correctness_of_biology_riddle(llm_response)
print(result)