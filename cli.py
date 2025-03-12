# import argparse
# import csv
# import sys
# import os
# from api import fetch_papers

# def main():
#     parser = argparse.ArgumentParser(description="Fetch research papers from PubMed.")
#     parser.add_argument("query", type=str, nargs="?", default="cancer treatment", help="Search query for PubMed")
#     parser.add_argument("-f", "--file", type=str, default="output.csv", help="Save results to a CSV file (default: output.csv)")
#     args = parser.parse_args()

#     print(f"Fetching papers from PubMed for query: {args.query}...")

#     # Ensure UTF-8 encoding for Windows output
#     sys.stdout.reconfigure(encoding="utf-8")

#     papers = fetch_papers(args.query)

#     file_path = os.path.abspath(args.file)  # Ensure full path is used
#     with open(file_path, "w", newline="", encoding="utf-8") as f:
#         writer = csv.DictWriter(f, fieldnames=["PubmedID", "Title", "Publication Date", "Non-academic Author(s)", "Company Affiliation(s)", "Corresponding Author Email"])
#         writer.writeheader()
#         writer.writerows(papers)
    
#     print(f"Results saved to {file_path}")

# if __name__ == "__main__":
#     main()

import argparse
import csv
import sys
from api import fetch_papers

def main():
    parser = argparse.ArgumentParser(description="Fetch research papers from PubMed.")
    parser.add_argument("query", type=str, help="Search query for PubMed")
    parser.add_argument("-f", "--file", type=str, help="Save results to a CSV file")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode")
    args = parser.parse_args()

    print(f"Fetching papers from PubMed for query: {args.query}...")

    # Set default encoding to UTF-8 (for Windows compatibility)
    sys.stdout.reconfigure(encoding="utf-8")

    # Fetch papers with debug mode if enabled
    papers = fetch_papers(args.query, debug=args.debug)

    if args.file:
        with open(args.file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["PubmedID", "Title", "Publication Date", 
                                                   "Non-academic Author(s)", "Company Affiliation(s)", 
                                                   "Corresponding Author Email"])
            writer.writeheader()
            writer.writerows(papers)
        print(f"Results saved to {args.file}")
    else:
        for paper in papers:
            print(paper)  # Print results to console

if __name__ == "__main__":
    main()
