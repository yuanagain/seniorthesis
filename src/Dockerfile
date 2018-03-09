FROM python:3


WORKDIR /usr/src/app

# Update
#RUN apk add --update python py-pip


COPY . .

# Install app dependencies
RUN python -m pip install --upgrade pip

RUN pip install --no-cache-dir -r requirements.txt

# Bundle app source
CMD ["python", "./experiment_mintrick3b.py"]

