import sys

def run_security_audit():
    """
    An interactive checklist to audit physical security controls for DNS infrastructure.
    This helps implement the Physical Security Guard Solution (F).
    """
    questions = {
        "Access Control": "Are data centers protected by multi-factor, role-based access controls? (y/n): ",
        "Visitor Management": "Is there a strict visitor management policy with escorted access and detailed logging? (y/n): ",
        "Surveillance": "Is the facility monitored 24/7 by high-resolution, tamper-proof CCTV? (y/n): ",
        "Intrusion Detection": "Are environmental and intrusion detection systems (e.g., motion sensors) in place? (y/n): ",
        "Hardware Security": "Are server racks secured with high-security locks and tamper-evident seals? (y/n): ",
        "Secure Deployment": "Are DNS servers deployed in secure, dedicated vaults or cages? (y/n): ",
        "Security Audits": "Are regular, independent physical security audits and penetration tests conducted? (y/n): ",
        "Hardware Lifecycle": "Is there a secure hardware lifecycle management process (from procurement to disposal)? (y/n): ",
        "Employee Vetting": "Do all personnel with physical access undergo comprehensive background checks? (y/n): ",
        "Incident Response": "Is there a documented and tested incident response plan for physical security breaches? (y/n): "
    }

    score = 0
    total_questions = len(questions)
    recommendations = []

    print("--- Physical Security Audit for DNS Infrastructure ---")
    print("Please answer the following questions with 'y' for yes or 'n' for no.\n")

    for control, question in questions.items():
        while True:
            answer = input(question).lower().strip()
            if answer in ['y', 'yes']:
                score += 1
                break
            elif answer in ['n', 'no']:
                recommendations.append(f"- Implement robust '{control}' measures.")
                break
            else:
                print("Invalid input. Please enter 'y' or 'n'.")

    security_level = (score / total_questions) * 100

    print("\n--- Audit Report ---")
    print(f"Audit Complete. You answered 'yes' to {score} out of {total_questions} questions.")
    print(f"Your physical security posture is at {security_level:.2f}%.")

    if security_level == 100:
        print("\nResult: Excellent! Your physical security controls appear to be comprehensive.")
    elif security_level >= 70:
        print("\nResult: Good. Your physical security is strong, but there are areas for improvement.")
    elif security_level >= 40:
        print("\nResult: Moderate. Significant gaps exist in your physical security controls.")
    else:
        print("\nResult: Critical. Your infrastructure is highly vulnerable to physical threats.")

    if recommendations:
        print("\nRecommendations for Improvement:")
        for rec in recommendations:
            print(rec)

if __name__ == '__main__':
    run_security_audit()