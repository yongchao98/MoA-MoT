import ipaddress

def calculate_summary_acl():
    """
    Calculates and prints the smallest summary ACL entry for a given list of networks.
    An ACL entry consists of a network address and a wildcard mask.
    """
    # Define the networks to be summarized
    network1 = ipaddress.ip_network("172.20.96.0/19")
    network2 = ipaddress.ip_network("172.20.128.0/19")

    # The ipaddress.collapse_addresses function finds the minimal representation
    # for a list of IP networks. For two adjacent, summarizable networks,
    # it will return a single supernet.
    summary_list = list(ipaddress.collapse_addresses([network1, network2]))

    # We expect a single summary network as the result
    if not summary_list or len(summary_list) > 1:
        print("Error: Could not summarize the given networks into a single entry.")
        return

    summary_network = summary_list[0]

    # The summary network address is the base IP of the supernet.
    summary_address = summary_network.network_address

    # The wildcard mask is the inverse of the subnet mask.
    # The ipaddress library provides this directly as the 'hostmask'.
    wildcard_mask = summary_network.hostmask

    # Print the final result, showing each number in the address and mask.
    # The format is: network_octet1 network_octet2 ... wildcard_octet1 ...
    addr_octets = str(summary_address).split('.')
    mask_octets = str(wildcard_mask).split('.')

    print(f"{addr_octets[0]}.{addr_octets[1]}.{addr_octets[2]}.{addr_octets[3]} {mask_octets[0]}.{mask_octets[1]}.{mask_octets[2]}.{mask_octets[3]}")

# Execute the function to find and print the ACL entry
calculate_summary_acl()