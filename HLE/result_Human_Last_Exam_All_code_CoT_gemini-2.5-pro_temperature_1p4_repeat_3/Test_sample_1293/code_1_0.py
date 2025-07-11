import dns.resolver
import threading

def run_verification():
    """
    This script demonstrates the Multi-Modal Verification Process for DNS security.
    It queries multiple independent DNS resolvers for the same domain name and
    compares their responses to detect inconsistencies, which could indicate an
    Adversary-in-the-Middle (AitM) attack.
    """
    # Step 1: Define multiple DNS resolvers to act as independent verification channels.
    resolvers_to_query = {
        'Google': '8.8.8.8',
        'Cloudflare': '1.1.1.1',
        'Quad9': '9.9.9.9',
        'OpenDNS': '208.67.222.222'
    }
    # Step 2: Define the target domain name to verify.
    domain_to_check = 'www.iana.org'
    results = {}
    lock = threading.Lock()

    def query_dns_resolver(name, ip):
        """Queries a specific DNS resolver and stores the result."""
        try:
            resolver = dns.resolver.Resolver(configure=False)
            resolver.nameservers = [ip]
            resolver.timeout = 3
            resolver.lifetime = 3
            answer = resolver.resolve(domain_to_check, 'A')
            # Store a sorted set of IPs for consistent comparison
            ips = set(sorted([rdata.to_text() for rdata in answer]))
            with lock:
                results[name] = ips
        except Exception as e:
            with lock:
                results[name] = f"Error: {e.__class__.__name__}"

    # Step 3: Use threads for simultaneous queries.
    threads = []
    print(f"[*] Starting Multi-Modal DNS Verification for domain: {domain_to_check}\n")

    for name, ip in resolvers_to_query.items():
        thread = threading.Thread(target=query_dns_resolver, args=(name, ip))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    # Step 4 & 5: Collect and report results from each channel.
    print("--- Individual Channel Verification Report ---")
    all_ip_sets = []
    has_errors = False
    for name, response in results.items():
        resolver_ip = resolvers_to_query[name]
        if isinstance(response, set):
            # For the output, we show the numbers (IPs) that lead to the conclusion
            print(f"[+] Result from {name} ({resolver_ip}): {sorted(list(response))}")
            all_ip_sets.append(response)
        else:
            print(f"[-] Result from {name} ({resolver_ip}): {response}")
            has_errors = True
    print("-" * 44)

    # Step 6: Analyze the collected data and provide a final conclusion.
    if has_errors:
        print("\n[!] FINAL CONCLUSION: Verification Inconclusive.")
        print("    Reason: One or more DNS queries failed. This could be a network issue or a targeted DoS on a resolver.")
    elif not all_ip_sets:
        print("\n[!] FINAL CONCLUSION: Verification Failed.")
        print("    Reason: No valid DNS responses were received from any channel.")
    else:
        # Check if all collected sets of IPs are identical
        first_set = all_ip_sets[0]
        is_consistent = all(s == first_set for s in all_ip_sets)

        if is_consistent:
            print("\n[SUCCESS] FINAL CONCLUSION: DNS Resolution is Consistent and Verified.")
            print(f"    All channels agree on the IP set for {domain_to_check}: {sorted(list(first_set))}")
        else:
            print("\n[CRITICAL ALERT] FINAL CONCLUSION: DNS Inconsistency Detected! Potential AitM Attack!")
            print("    Reason: Different resolvers returned different IP addresses for the same domain.")
            unique_sets = {}
            for i, s in enumerate(all_ip_sets):
                # Group resolvers by the IP set they returned
                s_tuple = tuple(sorted(list(s)))
                if s_tuple not in unique_sets:
                    unique_sets[s_tuple] = []
                # Find the name of the resolver that gave this set
                for name, res_set in results.items():
                    if res_set == s:
                        unique_sets[s_tuple].append(name)

            print("    Discrepant IP sets found:")
            for i, (ip_tuple, resolver_names) in enumerate(unique_sets.items()):
                print(f"      -> Set {i+1} from {', '.join(resolver_names)}: {list(ip_tuple)}")

if __name__ == '__main__':
    run_verification()