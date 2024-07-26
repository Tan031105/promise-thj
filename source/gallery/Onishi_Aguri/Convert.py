def urls_to_markdown():
    import sys
    import re
    
    print("请逐行输入URL，输入空行结束输入：")
    
    urls = []
    while True:
        line = sys.stdin.readline().strip()
        if line == "":
            break
        urls.append(line)
    
    markdown_links = []
    for i, url in enumerate(urls):
        match = re.match(r'^https?://', url)
        if match:
            markdown_links.append(f"![]({url})")
        else:
            print(f"无效的URL：{url}")
    
    for link in markdown_links:
        print(link)

if __name__ == "__main__":
    urls_to_markdown()
