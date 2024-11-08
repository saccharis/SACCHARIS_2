import platform
import re
import socket
import subprocess
import time
from urllib.parse import urlparse


def get_dns_servers():
    os_type = platform.system()

    try:
        if os_type in ["Linux", "Darwin"]:  # Darwin is macOS
            dns_servers = []
            with open('/etc/resolv.conf') as file:
                for line in file:
                    if line.startswith("nameserver"):
                        dns_servers.append(line.strip().split()[1])
            return dns_servers

        elif os_type == "Windows":
            output = subprocess.check_output("ipconfig /all", shell=True, text=True)

            # Regex to match DNS Servers and capture the subsequent IPs
            dns_matches = re.findall(r'(?:\d+\.\d+\.\d+\.\d+)+', re.search(r'DNS Servers[ .]*:(?:\s*[\d\.]*\n\s*)+', output).group())

            return dns_matches
        else:
            return "Unsupported operating system for DNS server lookup."

    except Exception as e:
        return f"Could not retrieve DNS servers: {e}"


def resolve_hostname(url):
    try:
        hostname = str(urlparse(url).hostname)
        if not hostname:
            return None, "Invalid URL; hostname could not be determined."

        start_time = time.time()
        ip_address = socket.gethostbyname(hostname)
        lookup_time = time.time() - start_time
        return ip_address, lookup_time
    except socket.gaierror as e:
        return None, f"DNS lookup failed: {e}"
