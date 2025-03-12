import ssl
import certifi
import json
from Bio import Entrez
from typing import List, Dict

# Force Python to use certifi's certificate
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

Entrez.email = "your_email@example.com"  # Replace with your email

def fetch_papers(query: str, max_results: int = 10, debug: bool = False) -> List[Dict]:
    """Fetches research papers from PubMed and extracts necessary details."""
    if debug:
        print(f"[DEBUG] Searching PubMed for: {query}")

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        paper_ids = record["IdList"]
        papers = []

        for paper_id in paper_ids:
            handle = Entrez.efetch(db="pubmed", id=paper_id, rettype="xml", retmode="xml")
            details = Entrez.read(handle)
            handle.close()

            if debug:
                print(f"\n=== [DEBUG] Raw API response for {paper_id} ===\n")
                print(json.dumps(details, indent=2))  # Pretty-print JSON
            
            try:
                article = details["PubmedArticle"][0]["MedlineCitation"]["Article"]
                authors_list = article.get("AuthorList", [])

                # Extract title
                title = article.get("ArticleTitle", "N/A")

                # Extract publication date
                pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
                year = pub_date.get("Year", "N/A")
                month = pub_date.get("Month", "")
                day = pub_date.get("Day", "")
                publication_date = f"{year}-{month}-{day}".strip("-")

                # Extract authors & affiliations
                non_academic_authors = []
                company_affiliations = []
                corresponding_email = "N/A"

                for author in authors_list:
                    last_name = author.get("LastName", "")
                    fore_name = author.get("ForeName", "")
                    full_name = f"{fore_name} {last_name}".strip()

                    # Extract affiliation
                    affiliations = author.get("AffiliationInfo", [])
                    for aff in affiliations:
                        aff_text = aff.get("Affiliation", "")
                        if any(keyword in aff_text.lower() for keyword in ["inc.", "ltd", "biotech", "pharmaceutical", "company"]):
                            company_affiliations.append(aff_text)
                        else:
                            non_academic_authors.append(full_name)

                    # Extract email (if present)
                    if "@" in aff_text:
                        corresponding_email = aff_text

                paper_info = {
                    "PubmedID": paper_id,
                    "Title": title,
                    "Publication Date": publication_date,
                    "Non-academic Author(s)": "; ".join(set(non_academic_authors)) or "N/A",
                    "Company Affiliation(s)": "; ".join(set(company_affiliations)) or "N/A",
                    "Corresponding Author Email": corresponding_email
                }

                papers.append(paper_info)

            except (KeyError, IndexError):
                if debug:
                    print(f"[DEBUG] Failed to parse data for PubmedID: {paper_id}")
                continue

        return papers

    except Exception as e:
        print(f"Error fetching papers: {e}")
        return []

