from flask import Flask, request, render_template
from align import Algorithm, Aligner
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

app = Flask(__name__)

@app.route('/')
def main():
    return homepage()


def homepage():
    html = 'Sequence <form method="POST"><input name="sequence">'
    html += '<select name="algorithm">'
    for algorithm in Algorithm:
        html += '<option value="{}]">{}</option>'.format(str(algorithm)[10:], str(algorithm)[10:])
    html += '</select>'
    html += '<input type="submit" value="Align"></form>'

    return html


@app.route('/', methods=['POST'])
def process():
    filename = './sequence.fasta'
    sequence = request.form['sequence']
    algorithm = request.form['algorithm']

    seq = SeqRecord(Seq(sequence), id="test")
    SeqIO.write(seq, filename, "fasta")
    aligner = Aligner('./db/db.txt')

    return homepage() + aligner.align(filename, Algorithm[algorithm[:-1]])


if __name__ == '__main__':
    app.run()
