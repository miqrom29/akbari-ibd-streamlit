FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 10000

HEALTHCHECK CMD curl --fail http://localhost:10000/_stcore/health

ENTRYPOINT ["streamlit", "run", "app.py", "--server.port=10000", "--server.address=0.0.0.0"]
