FROM python:3.7-bullseye
COPY requirments.txt /tmp
RUN pip install --no-cache-dir -r /tmp/requirments.txt && rm /tmp/requirments.txt
COPY . /app
WORKDIR /app
CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
