# 📚 PubMed Research Paper Fetcher

## 🚀 Overview

PubMed Research Paper Fetcher is a Python-based tool that allows users to fetch research papers from PubMed based on a given query. It extracts key information such as paper title, publication date, authors, affiliations, and corresponding author email, and saves the results in a CSV file.

## 🔥 Features

- Extracts key metadata such as **PubmedID, Title, Authors, Publication Date, Affiliations, and Email**.
- Saves results in a **CSV file** for easy access.

## 🛠 Installation

1. **Clone the Repository**
   ```sh
   git clone https://github.com/Saurajit-Nayak/project_pubmed_fetcher.git
   cd project_pubmed_fetcher
   ```
2. **Set Up Virtual Environment (Optional but Recommended)**
   ```sh
   python -m venv venv
   source venv/bin/activate  # On Windows use: venv\Scripts\activate
   ```
3. **Install Dependencies**
   ```sh
   pip install -r requirements.txt
   ```

## 🏃‍♂️ Usage

### Run from CLI

#### 1️⃣ Fetch and Print Papers to Console

```sh
python cli.py "cancer treatment"
```

#### 2️⃣ Save Results to a CSV File

```sh
python cli.py "cancer treatment" -f results.csv
```

#### 3️⃣ Enable Debug Mode

```sh
python cli.py "cancer treatment" -d
```

## 📝 Example Output

| PubmedID | Title                                                      | Publication Date | Authors    | Company Affiliation | Corresponding Email |
| -------- | ---------------------------------------------------------- | ---------------- | ---------- | ------------------- | ------------------- |
| 40066698 | Single-cell profiling reveals the intratumor heterogeneity | 2025-Mar-11      | Liang Weng | Key Lab Oncology    | N/A                 |

## 📌 Project Structure

```
project_pubmed_fetcher/
│── api.py              # Fetches PubMed papers
│── cli.py              # Command-line interface
│── requirements.txt    # Dependencies
│── README.md           # Project documentation
```

## 🤝 Contribution

1. **Create a new branch** (`git checkout -b feature-branch`)
2. **Make your changes**
3. **Commit your changes** (`git commit -m "Added new feature"`)
4. **Push to GitHub** (`git push origin feature-branch`)
5. **Open a Pull Request** 🚀

## 📬 Contact

For any questions or issues, feel free to reach out at **[saurajitnayak095@gmail.com](mailto\:saurajitnayak095@gmail.com)**.

---

🚀 **Happy Researching!**

