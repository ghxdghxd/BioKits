#!/bin/bash
SSL_PATH="/etc/gitlab/ssl"
name="gitlab"

#建立认证目录
mkdir $SSL_PATH
chmod 700 $SSL_PATH

#建立证书 自签名证书
#创建 Private Key
openssl genrsa -des3 -out $SSL_PATH/$name\.key 2048
openssl req -new -key $SSL_PATH/$name.key -out $SSL_PATH/$name\.csr

# openssl req -nodes -newkey rsa:2048 -keyout $SSL_PATH/$name\.key -out $SSL_PATH/$name\.csr

# Enter Country Name US
# Enter State or Province Full Name
# Enter City Name
# Enter Organization Name
# Enter Company Name
# Enter Organizational Unit Name
# Enter server hostname i.e. URL gitlab.domain.com
# Enter Admin Email Address
# Skip Challenge Password (Hit Enter)
# Skip Optional Company Name (Hit Enter)

#移除Private Key 中的密码短语
cp -v $SSL_PATH/$name\.{key,original}
openssl rsa -in $SSL_PATH/$name\.original -out $SSL_PATH/$name\.key
rm -v $SSL_PATH/$name\.original
#创建证书
openssl x509 -req -days 1460 -in $SSL_PATH/$name\.csr -signkey $SSL_PATH/$name\.key -out $SSL_PATH/$name\.crt
# 移除证书请求文件
rm -v $SSL_PATH/$name\.csr
#设置文件权限
chmod 600 $SSL_PATH/$name\.*