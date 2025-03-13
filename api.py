import ssl
import certifi
import json
import time
from Bio import Entrez
from typing import List, Dict

ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())


Entrez.email = "saurajitnayak095@gmail.com"

def fetch_papers(query: str, max_results: int = 5, debug: bool = False) -> List[Dict]:
    """Fetches research papers from PubMed and extracts necessary details."""
    if debug:
        print(f"[DEBUG] Searching PubMed for: {query}")

    try:
       
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        paper_ids = record.get("IdList", [])
        if debug:
            print(f"[DEBUG] Total results from PubMed: {len(paper_ids)}")

        papers = []

        for paper_id in paper_ids:
            time.sleep(1)  # Prevents API rate-limiting

            # Step 2: Fetch paper details
            handle = Entrez.efetch(db="pubmed", id=paper_id, rettype="xml", retmode="xml")
            details = Entrez.read(handle)
            handle.close()

            if debug:
                print(f"\n=== [DEBUG] Raw API response for {paper_id} ===\n")
                print(json.dumps(details, indent=2))  

            try:
                article = details["PubmedArticle"][0]["MedlineCitation"]["Article"]
                authors_list = article.get("AuthorList", [])

               
                title = article.get("ArticleTitle", "N/A")

               
                pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
                year = pub_date.get("Year", "Unknown")
                month = pub_date.get("Month", "Unknown")
                day = pub_date.get("Day", "Unknown")
                publication_date = f"{year}-{month}-{day}".replace("Unknown-", "")

                
                non_academic_authors = []
                company_affiliations = []
                corresponding_email = "N/A"

                for author in authors_list:
                    last_name = author.get("LastName", "")
                    fore_name = author.get("ForeName", "")
                    full_name = f"{fore_name} {last_name}".strip()

                    affiliations = author.get("AffiliationInfo", [])
                    for aff in affiliations:
                        aff_text = aff.get("Affiliation", "")

                       
                        if any(keyword in aff_text.lower() for keyword in ["inc.", "ltd", "biotech", "pharmaceutical", "company"]):
                            company_affiliations.append(aff_text)
                        else:
                            non_academic_authors.append(full_name)

                    
                    email = author.get("ElectronicAddress", "")
                    if "@" in email:
                        corresponding_email = email

                
                if not non_academic_authors:
                    non_academic_authors.append("N/A")

             
                paper_info = {
                    "PubmedID": paper_id,
                    "Title": title,
                    "Publication Date": publication_date,
                    "Non-academic Author(s)": "; ".join(set(non_academic_authors)),
                    "Company Affiliation(s)": "; ".join(set(company_affiliations)) or "N/A",
                    "Corresponding Author Email": corresponding_email
                }

                papers.append(paper_info)

            except (KeyError, IndexError) as e:
                if debug:
                    print(f"[DEBUG] Failed to parse data for PubmedID {paper_id}: {e}")
                continue

        return papers

    except Exception as e:
        print(f"Error fetching papers: {e}")
        return []


if __name__ == "__main__":
    query = "Melatonin cancer"
    papers = fetch_papers(query, max_results=10, debug=True)
    for paper in papers:
        print(paper)


