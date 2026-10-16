FROM python:3.12-slim

# Hugging Face security requirement
RUN useradd -m -u 1000 user
USER user
ENV PATH="/home/user/.local/bin:$PATH"

WORKDIR /app

# Copy the entire Oligo_Designer_Tool folder into the container
COPY --chown=user . /app

# Install the dependencies
RUN pip install --no-cache-dir -r requirements.txt -r requirements-server.txt

EXPOSE 7860

# Run the FastAPI server directly from this folder
CMD ["uvicorn", "fastapi_backend:app", "--host", "0.0.0.0", "--port", "7860"]
