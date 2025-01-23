```bash
python -m venv venv
source ./venv/bin/activate
pip install .
```

proteins 主要负责蛋白质相关的 CRUD 操作。

stats 主要负责统计相关的 CRUD 操作。

```bash
python manage.py makemigrations
python manage.py migrate
python manage.py load_proteins -p /path/to/data.csv
```
