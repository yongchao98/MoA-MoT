import ipaddress

def calculate_summary_acl(networks):
    """
    Calculates the smallest single ACL entry (supernet) to cover a list of networks.

    Args:
        networks (list): A list of network strings in CIDR format.

    Returns:
        tuple: A tuple containing the summary network address and wildcard mask as strings.
    """
    if not networks:
        return None, None

    # Convert string networks to ipaddress objects
    ip_networks = [ipaddress.ip_network(n, strict=False) for n in networks]

    # Find the absolute first and last IP addresses in the combined range
    first_ip = min(net.network_address for net in ip_networks)
    last_ip = max(net.broadcast_address for net in ip_networks)

    # The easiest way to get the summary is to use ipaddress.summarize_address_range
    # This returns a generator of the smallest list of networks covering the range.
    # For this specific problem, it will yield one summary network.
    summary_net = list(ipaddress.summarize_address_range(first_ip, last_ip))[0]

    summary_address = summary_net.network_address
    wildcard_mask = summary_net.hostmask

    return str(summary_address), str(wildcard_mask)

# The networks given in the problem
network_list = ['172.20.96.0/19', '172.20.128.0/19']

# Calculate the summary ACL entry
summary_addr, wildcard = calculate_summary_acl(network_list)

# Print the final result
if summary_addr and wildcard:
    print(f"The smallest appropriate ACL entry is:")
    # Output each number in the final equation as requested
    print(f"{summary_addr} {wildcard}")
else:
    print("Could not calculate summary route.")
