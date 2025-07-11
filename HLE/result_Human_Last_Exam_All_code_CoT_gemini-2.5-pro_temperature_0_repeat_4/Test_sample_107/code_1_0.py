import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest ACL entry (supernet) for a given list of networks.
    """
    # The input networks from the problem description
    network_strings = [
        '172.20.96.0/19',
        '172.20.128.0/19'
    ]

    # Convert the string representations into ipaddress network objects
    networks = [ipaddress.ip_network(net) for net in network_strings]

    # The collapse_addresses function summarizes a list of networks
    # into the smallest possible list of covering networks.
    # For these two networks, it will return a single summary network (supernet).
    summary_list = list(ipaddress.collapse_addresses(networks))

    # The result is the first (and only) network in the summarized list
    supernet = summary_list[0]

    # The network address of the supernet
    summary_address = supernet.network_address

    # The wildcard mask is the inverse of the netmask, which the ipaddress
    # library provides as the 'hostmask' attribute.
    wildcard_mask = supernet.hostmask

    # Print the final ACL entry in the format: network_address wildcard_mask
    print(f"{summary_address} {wildcard_mask}")

if __name__ == "__main__":
    calculate_summary_acl()