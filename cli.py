
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

    sys.stdout.reconfigure(encoding="utf-8")

    papers = fetch_papers(args.query, debug=args.debug)

    if not papers:
        print("No relevant papers found. Try a different query.")
        return

    if args.debug:
        print(f"[DEBUG] Fetched {len(papers)} papers from PubMed.")

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
            print(paper)

if __name__ == "__main__":
    main()
