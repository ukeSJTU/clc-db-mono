这个项目采用了 monorepo 的结构，前端 nextjs 后端 django。

后端配置

先配置 python 环境

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python manage.py migrate
```

如果是第一次运行，需要导入数据库数据。

```bash
python manage.py makemigrations
python manage.py migrate
python manage.py load_proteins -p /path/to/data.csv
```

然后是根据 sdf 生成 2D 图片。

```bash
python manage.py makemigrations
python manage.py migrate
python manage.py load_proteins -p /path/to/data.csv
```

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python manage.py migrate
guvicorn clcdb.asgi:application --bind
```

下面是项目的部署方式。先拷贝 nginx 配置文件。

```bash
sudo cp ./clcdb.conf /etc/nginx/sites-available/
sudo ln -s /etc/nginx/sites-available/clcdb.conf /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```
