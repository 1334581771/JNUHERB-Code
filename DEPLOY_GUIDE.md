# JNUHERB 京东云轻量云主机完整部署教程（小白版）

> 适用环境：京东云轻量云主机 2核4G / Ubuntu 22.04+ / 已安装 nginx
> 目标：把 JNUHERB 项目部署成可通过浏览器访问的网页，并支持本地修改后快速更新
> 部署方式：Docker（推荐）
> 你的服务器信息：IP `111.228.21.145`，用户名 `root`

---

## 一、准备工作

### 1.1 你需要准备的信息

- [ ] 云主机公网 IP：`111.228.21.145`
- [ ] root 登录密码
- [ ] 本地项目路径（例如 Windows 上 `D:\herbtar\jnuherb`）
- [ ] 域名（可选，没有也能用 IP 访问）
- [ ] SSH 客户端：
  - Windows：推荐 [Termius](https://termius.com/)、[FinalShell](https://www.hostbuf.com/) 或 PowerShell
  - Mac/Linux：系统自带“终端”

### 1.2 用 SSH 登录服务器

```bash
ssh root@111.228.21.145
```

输入密码后看到 `root@xxx:~#` 即表示登录成功。

> 如果登录失败，检查：IP、用户名、密码、京东云防火墙是否放行 22 端口。

---

## 二、Docker 部署（推荐）

### 2.1 安装 Docker

京东云访问 Docker 官方海外源可能失败，因此**推荐用国内源安装**。

```bash
# 更新系统
sudo apt update && sudo apt upgrade -y

# 安装依赖
sudo apt install -y ca-certificates curl gnupg lsb-release

# 添加阿里云 Docker GPG 密钥
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://mirrors.aliyun.com/docker-ce/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg

# 添加阿里云 Docker 软件源
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://mirrors.aliyun.com/docker-ce/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# 安装 Docker
sudo apt update
sudo apt install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin

# 验证
sudo docker --version
```

> 如果国内源也失败，可改用 Ubuntu 自带源：`sudo apt install -y docker.io docker-compose`

### 2.2 验证 docker compose

新版 Docker 的命令是 `docker compose`（中间有空格），不是旧的 `docker-compose`。

```bash
sudo docker compose version
```

看到版本号即可。

### 2.3 配置 Docker 国内镜像加速

因为 Docker Hub 在国外，京东云直接拉取容易失败，必须配置国内镜像加速：

```bash
sudo mkdir -p /etc/docker
sudo tee /etc/docker/daemon.json <<-'EOF'
{
  "registry-mirrors": [
    "https://docker.mirrors.ustc.edu.cn",
    "https://hub-mirror.c.163.com",
    "https://mirror.baidubce.com",
    "https://docker.m.daocloud.io"
  ]
}
EOF

sudo systemctl daemon-reload
sudo systemctl restart docker
```

验证是否生效：

```bash
sudo docker info | grep -A 10 "Registry Mirrors"
```

### 2.4 上传 JNUHERB 项目到服务器

#### 方式一：rsync（推荐，只传改动）

在本地 PowerShell / 终端执行：

```powershell
rsync -avz --delete `
  --exclude='src/jnuherb/data/*.csv' `
  --exclude='src/jnuherb/data/jnuherb.db' `
  --exclude='.git' `
  --exclude='__pycache__' `
  --exclude='*.pyc' `
  /d/herbtar/jnuherb/ `
  root@111.228.21.145:/opt/jnuherb/
```

> 如果 Windows 没有 rsync，可用 WinSCP 或安装 Git for Windows / WSL。

#### 方式二：WinSCP（图形化，小白友好）

1. 下载 [WinSCP](https://winscp.net/)
2. 协议选 **SFTP**，主机 `111.228.21.145`，用户名 `root`
3. 左边本地打开 `D:\herbtar\jnuherb`，右边服务器打开 `/opt/jnuherb`
4. 把文件拖过去覆盖

#### 方式三：scp

```bash
scp -r /d/herbtar/jnuherb root@111.228.21.145:/opt/
```

### 2.5 启动 JNUHERB

```bash
cd /opt/jnuherb
sudo docker compose down
sudo docker compose up -d --build
```

第一次构建会比较慢（10–90 分钟），因为需要下载 Python 基础镜像和安装依赖。

### 2.6 检查启动状态

```bash
sudo docker compose ps
sudo docker compose logs -f
```

正常应看到：

```text
NAME      STATUS          PORTS
jnuherb   Up              0.0.0.0:8501->8501/tcp
```

日志中应出现：

```text
App URL: http://0.0.0.0:8501
```

### 2.7 放行防火墙端口

登录京东云控制台 → 找到轻量云主机 → 防火墙/安全组 → 添加规则：

| 协议 | 端口 | 来源 |
|---|---|---|
| TCP | 8501 | 0.0.0.0/0 |
| TCP | 80 | 0.0.0.0/0 |
| TCP | 443 | 0.0.0.0/0 |

如果只打算用 `IP:8501` 访问，先只放行 8501 也行。

### 2.8 浏览器访问

```text
http://111.228.21.145:8501
```

看到 JNUHERB 界面即部署成功。

---

## 三、用 nginx 配置域名/IP 访问（可选但推荐）

### 3.1 域名解析（可选）

把域名 A 记录指向 `111.228.21.145`，等待 5–30 分钟生效。

### 3.2 配置 nginx

```bash
sudo nano /etc/nginx/sites-available/jnuherb
```

粘贴以下内容（用 IP 访问就写 IP，有域名就写域名）：

```nginx
server {
    listen 80;
    server_name 111.228.21.145;

    location / {
        proxy_pass http://127.0.0.1:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 420s;
        proxy_connect_timeout 420s;
        proxy_send_timeout 420s;
    }
}
```

启用：

```bash
sudo ln -s /etc/nginx/sites-available/jnuherb /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

然后访问 `http://111.228.21.145` 即可（不用带端口号）。

---

## 四、配置 HTTPS（有域名时推荐）

```bash
sudo apt install -y certbot python3-certbot-nginx
sudo certbot --nginx -d 你的域名.com
```

按提示输入邮箱、同意条款。成功后自动启用 HTTPS。

访问：`https://你的域名.com`

---

## 五、本地修改后上传到服务器并更新

这是你的长期工作流：本地改代码 → 上传服务器 → 重启 Docker。

### 5.1 上传改动

```powershell
rsync -avz --delete `
  --exclude='src/jnuherb/data/*.csv' `
  --exclude='src/jnuherb/data/jnuherb.db' `
  --exclude='.git' `
  --exclude='__pycache__' `
  --exclude='*.pyc' `
  /d/herbtar/jnuherb/ `
  root@111.228.21.145:/opt/jnuherb/
```

> 只改了一个文件也可以用 scp：
> ```bash
> scp /d/herbtar/jnuherb/src/jnuherb/data/streamlit_app.py root@111.228.21.145:/opt/jnuherb/src/jnuherb/data/
> ```

### 5.2 登录服务器重启

```bash
ssh root@111.228.21.145
cd /opt/jnuherb
sudo docker compose down
sudo docker compose up -d --build
```

### 5.3 确认更新真的生效了（关键！）

很多新手会遇到一个问题：**本地代码改了，上传了，也重启了，但网页还是旧的样子**。这是因为 Docker 镜像有缓存，或者旧镜像没被覆盖。

#### 第一步：确认宿主机上的源码已更新

```bash
grep -n "boxplot" /opt/jnuherb/src/jnuherb/data/streamlit_app.py
```

应该看到：

```python
ax.boxplot(data_groups)
ax.set_xticklabels(list(attrs))
```

#### 第二步：确认容器内的源码也已更新

```bash
sudo docker exec -it jnuherb grep -n "boxplot" /home/jnuherb/.local/lib/python3.11/site-packages/jnuherb/data/streamlit_app.py
```

如果这里还是旧的，说明 Docker 镜像没更新。

#### 第三步：如果容器内还是旧的，强制清理重建

```bash
cd /opt/jnuherb
sudo docker compose down
sudo docker image rm jnuherb:latest
sudo docker compose build --no-cache
sudo docker compose up -d
```

#### 第四步：应急方案——直接在容器内修改

如果临时要修一个小 bug，不想重新构建镜像，可以直接进容器改：

```bash
# 进入容器
sudo docker exec -it jnuherb bash

# 用 sed 替换有问题的代码（以 boxplot labels 参数为例）
sed -i 's/ax.boxplot(data_groups, labels=list(attrs))/ax.boxplot(data_groups)\n                ax.set_xticklabels(list(attrs))/' /home/jnuherb/.local/lib/python3.11/site-packages/jnuherb/data/streamlit_app.py

# 退出容器
exit

# 重启容器生效
sudo docker compose restart
```

> ⚠️ 注意：容器内修改只是临时的，容器重建后会丢失。最终还是要改本地源码并重新上传构建。

### 5.4 浏览器验证

```text
http://111.228.21.145:8501
```

运行一次分析，进入对应功能页面，确认修改已生效。

---

## 六、DeepSeek API Key 配置（可选）

AI 辨证功能需要 DeepSeek API Key。当前项目配置为在网页侧边栏手动填写，最灵活也最安全。

如果你希望服务器自动注入，可修改 `docker-compose.yml`，加入环境变量：

```yaml
environment:
  - DEEPSEEK_API_KEY=${DEEPSEEK_API_KEY:-}
env_file:
  - .env
```

然后在服务器创建 `.env`：

```bash
echo 'DEEPSEEK_API_KEY=sk-你的key' | sudo tee /opt/jnuherb/.env
```

> 注意：`.env` 文件已加入 `.gitignore`，不会提交到 git。

---

## 七、设置开机自动启动

配置 systemd 服务后，服务器重启会自动启动 JNUHERB，无需手动运行 Docker。

### 7.1 确保 Docker 本身开机自启

```bash
sudo systemctl enable docker
```

### 7.2 创建 JNUHERB 自启服务

```bash
sudo nano /etc/systemd/system/jnuherb.service
```

粘贴：

```ini
[Unit]
Description=JNUHERB Streamlit Application
Requires=docker.service
After=docker.service

[Service]
Type=oneshot
RemainAfterExit=yes
WorkingDirectory=/opt/jnuherb
# Use /bin/sh -c because systemd ExecStart does not handle "docker compose" subcommand reliably
ExecStart=/bin/sh -c 'docker compose up -d'
ExecStop=/bin/sh -c 'docker compose down'
TimeoutStartSec=0

[Install]
WantedBy=multi-user.target
```

### 7.3 启用并启动服务

```bash
sudo systemctl daemon-reload
sudo systemctl enable jnuherb
sudo systemctl start jnuherb
```

### 7.4 验证

```bash
sudo systemctl status jnuherb
sudo docker compose ps
```

### 7.5 常用管理命令

```bash
# 查看状态
sudo systemctl status jnuherb

# 手动重启
sudo systemctl restart jnuherb

# 手动停止
sudo systemctl stop jnuherb

# 查看日志
sudo journalctl -u jnuherb -f
```

---

## 八、常见问题排查

### 8.1 Docker 安装时报 `Could not handshake: Error in the pull function`

原因：连不上 Docker 官方海外源。  
解决：用本教程 2.1 节的国内源安装。

### 8.2 `docker-compose: command not found`

原因：新版 Docker 命令是 `docker compose`（空格）。  
解决：用 `sudo docker compose up -d`。

### 8.3 拉镜像时报 `connection refused` / `403 Forbidden`

原因：连不上 Docker Hub。  
解决：按 2.3 节配置国内镜像加速。

### 8.4 容器一直 `Restarting`

查看日志：

```bash
sudo docker compose logs -f
```

常见原因：

- **`--server.headless` 参数错误**：新版 Streamlit 需要 `--server.headless=true`
- **代码语法错误**：上传前先在本地跑通
- **数据文件缺失或路径错误**：检查 CSV 是否在 `src/jnuherb/data/`

### 8.5 网页打不开

依次检查：

```bash
sudo docker compose ps                    # 容器是否 Up
sudo ss -tlnp | grep 8501                # 8501 端口是否在监听
curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:8501  # 本地是否通
sudo ufw status                          # Ubuntu 防火墙是否放行
```

然后检查京东云控制台防火墙是否放行 8501/80/443。

### 8.6 报错 `Axes.boxplot() got an unexpected keyword argument 'labels'`

原因：matplotlib 新版（3.10+）把 `labels` 参数改名了。  
解决：把 `streamlit_app.py` 中：

```python
ax.boxplot(data_groups, labels=list(attrs))
```

改成：

```python
ax.boxplot(data_groups)
ax.set_xticklabels(list(attrs))
```

然后重新上传、重启 Docker。

### 8.7 上传后网页没有变化

这是更新代码时最常见的问题。根本原因是：**Docker 容器里运行的还是旧镜像**。

#### 排查步骤

1. **确认宿主机源码已更新**
   ```bash
   grep -n "boxplot" /opt/jnuherb/src/jnuherb/data/streamlit_app.py
   ```

2. **确认容器内源码也已更新**
   ```bash
   sudo docker exec -it jnuherb grep -n "boxplot" /home/jnuherb/.local/lib/python3.11/site-packages/jnuherb/data/streamlit_app.py
   ```

3. **如果容器内还是旧的，强制清理重建**
   ```bash
   cd /opt/jnuherb
   sudo docker compose down
   sudo docker image rm jnuherb:latest
   sudo docker compose build --no-cache
   sudo docker compose up -d
   ```

#### 实际案例：matplotlib `labels` 参数报错

你修改了本地 `streamlit_app.py`，上传后网页仍然报错：

```text
TypeError: Axes.boxplot() got an unexpected keyword argument 'labels'
```

排查发现宿主机上代码已改，但容器内还是旧的 `ax.boxplot(data_groups, labels=list(attrs))`。

解决方案就是上面的**强制清理重建**。重建后容器内代码变成：

```python
ax.boxplot(data_groups)
ax.set_xticklabels(list(attrs))
```

网页就不再报错了。

#### 为什么会出现这种情况

Docker 镜像由多层构成，build 时如果某一层（如 `COPY src/`）被判定为没有变化，就可能用缓存。加上 `--no-cache` 并删除旧镜像，可以确保 100% 重新构建。

#### 应急方案

如果暂时不想重建镜像，可以直接进入运行中的容器修改文件：

```bash
sudo docker exec -it jnuherb bash
sed -i 's/ax.boxplot(data_groups, labels=list(attrs))/ax.boxplot(data_groups)\n                ax.set_xticklabels(list(attrs))/' /home/jnuherb/.local/lib/python3.11/site-packages/jnuherb/data/streamlit_app.py
exit
sudo docker compose restart
```

> ⚠️ 注意：容器内修改是临时的，下次重建镜像会丢失。最终还是要改本地源码并重新上传构建。

---

## 九、命令速查表

| 场景 | 命令 |
|---|---|
| 登录服务器 | `ssh root@111.228.21.145` |
| 上传整个项目 | `rsync -avz --delete --exclude='src/jnuherb/data/*.csv' --exclude='src/jnuherb/data/jnuherb.db' /d/herbtar/jnuherb/ root@111.228.21.145:/opt/jnuherb/` |
| 上传单个文件 | `scp /d/herbtar/jnuherb/src/jnuherb/data/streamlit_app.py root@111.228.21.145:/opt/jnuherb/src/jnuherb/data/` |
| 停止容器 | `sudo docker compose down` |
| 启动/重建容器 | `sudo docker compose up -d --build` |
| 重启容器 | `sudo docker compose restart` |
| 查看状态 | `sudo docker compose ps` |
| 查看日志 | `sudo docker compose logs -f` |
| 测试本地 8501 | `curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:8501` |
| 重启 nginx | `sudo systemctl restart nginx` |
| 查看 nginx 配置 | `sudo nginx -t` |
| 查看系统防火墙 | `sudo ufw status` |

---

## 小服务器部署注意事项（2C4G 及以下）

如果你使用的是 2 核 4G 的轻量云主机，建议采取以下措施防止单个分析任务拖垮整个服务。
默认配置已允许最多 **2 个分析任务同时运行**，请根据实际并发量调整资源限制。

### 1. Docker 资源限制

`docker-compose.yml` 中已经配置了资源限制：

```yaml
deploy:
  resources:
    limits:
      cpus: '1.5'
      memory: 3G
    reservations:
      cpus: '0.5'
      memory: 1G
```

这确保 JNUHERB 容器最多使用 1.5 核 CPU 和 3G 内存，不会占满整个服务器。
**注意**：当 5 个分析任务同时运行时，1.5 核 / 3G 资源可能不足，建议升级到 4 核 8G
以上，或在 `docker-compose.yml` 中调高 `cpus` 和 `memory`。

### 2. 分析任务并发锁

JNUHERB 已内置文件锁机制：同一时刻最多允许 **2 个分析任务同时运行**。
如果服务器当前已有 2 个分析任务在运行，新用户点击 **Run Analysis** 会看到：

> "Server is busy with 2 analyses. Please wait a moment and try again."

如需修改并发数，请编辑 `src/jnuherb/data/streamlit_app.py` 中的常量：

```python
ANALYSIS_MAX_CONCURRENT = 2
```

改大该值后，务必同步增加服务器的 CPU / 内存配额，并调整 `docker-compose.yml`
中的资源限制，否则容易出现超时或 OOM。

### 3. 超时保护

单次分析最多运行 360 秒，超过后会自动停止并提示用户减少输入规模或调低 score threshold。这样即使用户关闭浏览器，服务器也不会被无限占用。

### 4. 输入规模限制

- Reverse 模式：最多 30 个目标蛋白
- Forward 模式：最多 30 个草药/方剂/单体
- Formula 名称模糊匹配：解析后最多 30 个方剂
- Forward 蛋白过滤：最多 50 个蛋白

### 5. nginx 超时设置

nginx 代理超时已设置为 300 秒，避免长连接一直挂起：

```nginx
proxy_read_timeout 300s;
proxy_connect_timeout 300s;
proxy_send_timeout 300s;
```

### 6. 应急处理

如果网页完全卡死、无法进入，登录服务器强制重启：

```bash
ssh root@111.228.21.145
cd /opt/jnuherb
sudo docker compose down
sudo docker compose up -d
```

如果怀疑有僵尸进程，可以：

```bash
sudo docker stats jnuherb          # 查看容器资源占用
sudo docker compose logs -f        # 查看实时日志
```

### 7. 分子对接（Docking）功能注意事项

网站的「🧪 Docking」栏目使用 AutoDock Vina 做单分子对接（exhaustiveness 8），
Docker 镜像已内置全部依赖（rdkit / meeko / vina，vina 在构建阶段编译）。

- **单次对接耗时**：2 核服务器上约 1–5 分钟，受体越大越慢
- **蛋白结构来源**：自动模式从 AlphaFold DB 在线下载（首次每个蛋白约几秒），
  下载结果在容器生命周期内缓存；容器重建后需重新下载
- **Ensembl → UniProt 映射**：打包在镜像内（`ensembl_uniprot.csv`），
  如数据库蛋白表更新，需重新运行 `scripts/build_uniprot_mapping.py` 并重建镜像
- **并发建议**：对接期间容器 CPU 占用较高，2C4G 服务器建议同时只进行
  1 个对接任务（用户连续点击会被 Streamlit 排队，一般不会冲突）

---

## 十、部署检查清单

- [ ] SSH 能登录服务器
- [ ] Docker 安装成功：`sudo docker --version`
- [ ] `docker compose` 可用：`sudo docker compose version`
- [ ] 国内镜像加速已配置
- [ ] 项目上传到 `/opt/jnuherb`
- [ ] `sudo docker compose up -d --build` 启动成功
- [ ] 京东云防火墙放行 8501（以及 80/443）
- [ ] 浏览器能访问 `http://111.228.21.145:8501`
- [ ] nginx 配置完成并能用 `http://111.228.21.145` 访问
- [ ] HTTPS 证书配置完成（可选）
- [ ] 开机自启服务配置完成
- [ ] DeepSeek API Key 已配置（如需 AI 功能）

---

## 十一、总结

对于小白来说，推荐的标准流程：

1. **首次部署**：安装 Docker → 配置镜像加速 → 上传项目 → `docker compose up -d --build`
2. **日常更新**：本地改代码 → rsync 上传 → 服务器上 `docker compose up -d --build`
3. **访问**：`http://111.228.21.145:8501`
4. **出问题**：先看 `sudo docker compose logs -f`

遇到问题把报错日志贴出来，就能快速定位。祝你部署顺利！🚀
